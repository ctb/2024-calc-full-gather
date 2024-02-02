#! /usr/bin/env python
# copied from https://github.com/ctb/2024-calc-full-gather commit 87d7d9ec8c67c
# for now.
"""

CTB TODO:
- deal with abundance stuff
- deal with threshold bp
- support the usual gather output? matches, matched hashes, etc.
- multiple databases...
"""
import argparse
import sys
import sourmash
import csv

from sourmash.search import GatherResult, format_bp
from sourmash.logging import print_results, set_quiet


def get_ident(name):
    return name.split(' ')[0]


def zipfile_load_ss_from_row(db, row):
    data = db.storage.load(row['internal_location'])
    sigs = sourmash.signature.load_signatures(data)

    return_sig = None
    for ss in sigs:
        if ss.md5sum() == row['md5']:
            assert return_sig is None # there can only be one!
            return_sig = ss

    if return_sig is None:
        raise ValueError("no match to requested row in db")
    return return_sig


def main():
    p = argparse.ArgumentParser()
    p.add_argument('query', help='query sketch')
    p.add_argument('database', help='database of sketches (zip file)')
    p.add_argument('fastgather_csv', help='output from fastgather')
    p.add_argument('--scaled', type=int, default=1000,
                   help='scaled value for comparison')
    p.add_argument('--threshold-bp', type=int, default=50000,
                   help='threshold for matches')
    p.add_argument('-o', '--output', default=None,
                   help='CSV output')
    p.add_argument('-q', '--quiet', default=False, action='store_true',
                   help='suppress output')
    p.add_argument('--estimate-ani-ci', default=False, action='store_true',
                   help='estimate ANI confidence intervals (default: False)')
    args = p.parse_args()

    set_quiet(args.quiet)

    db = sourmash.load_file_as_index(args.database)

    fastgather_results = []
    with open(args.fastgather_csv, 'r', newline='') as fp:
        r = csv.DictReader(fp)
        fastgather_results.extend(r)

    print(f"loaded {len(fastgather_results)} results.")
    for header in 'query_filename,rank,query_name,query_md5,match_name,match_md5,intersect_bp'.split(','):
        assert header in fastgather_results[0].keys()

    # find manifest entries => load directly? do we have that API?

    # or, do picklist?
    pl = sourmash.picklist.SignaturePicklist('prefetch',
                                             pickfile=args.fastgather_csv)
    _, dup_vals = pl.load()

    mf = db.manifest
    mf = mf.select_to_manifest(picklist=pl)

    # order rows by rank/order in gather result
    ident_to_row = {}
    for row in mf.rows:
        name = row['name']
        ident = get_ident(name)
        ident_to_row[ident] = row

    ordered_rows = []
    for n, gather_result in enumerate(fastgather_results):
        assert n == int(gather_result['rank'])
        ident = get_ident(gather_result['match_name'])
        mf_row = ident_to_row[ident]
        ordered_rows.append(mf_row)

    # guess ksize, get scaled - from first match
    first_ss = zipfile_load_ss_from_row(db, ordered_rows[0])
    ksize = first_ss.minhash.ksize
    scaled = max(args.scaled, first_ss.minhash.scaled)

    print(f"ksize={ksize}, scaled={scaled}")

    query_ss = sourmash.load_file_as_index(args.query)
    query_ss = query_ss.select(ksize=ksize, scaled=scaled)
    query_ss = list(query_ss.signatures())
    assert len(query_ss) == 1, query_ss
    query_ss = query_ss[0]
    assert query_ss.minhash.track_abundance, \
        "Query signatures must have abundance (for now)."

    orig_query_mh = query_ss.minhash.downsample(scaled=scaled)
    query_mh = orig_query_mh.to_mutable()

    orig_query_abunds = query_mh.hashes
    sum_abunds = sum(orig_query_abunds.values())

    # initialize output
    csv_writer = sourmash.sourmash_args.FileOutputCSV(args.output)
    outfp = csv_writer.open()
    result_writer = None

    # iterate over results, row by row
    screen_width = 80
    is_abundance = True
    sum_f_uniq_found = 0.
    found = False
    for rank, mf_row in enumerate(ordered_rows):
        best_match = zipfile_load_ss_from_row(db, mf_row)

        found_mh = best_match.minhash.downsample(scaled=scaled).flatten()
    
        n_weighted_missed = sum(( orig_query_abunds[k] for k in query_mh.hashes ))
        sum_weighted_found = sum_abunds - n_weighted_missed

        result = GatherResult(query_ss,
                              best_match,
                              cmp_scaled=scaled,
                              filename=args.database,
                              gather_result_rank=rank,
                              gather_querymh=query_mh,
                              ignore_abundance=False,
                              threshold_bp=args.threshold_bp,
                              orig_query_len=len(orig_query_mh),
                              orig_query_abunds=orig_query_abunds,
                              estimate_ani_ci=args.estimate_ani_ci,
                              sum_weighted_found=sum_weighted_found,
                              total_weighted_hashes=sum_abunds)

        query_mh.remove_many(found_mh)

        sum_f_uniq_found += result.f_unique_to_query

        if not found:                # first result? print header.
            if is_abundance:
                print_results("")
                print_results("overlap     p_query p_match avg_abund")
                print_results("---------   ------- ------- ---------")
            else:
                print_results("")
                print_results("overlap     p_query p_match")
                print_results("---------   ------- -------")

            found = True


        # print interim result & save in `found` list for later use
        pct_query = '{:.1f}%'.format(result.f_unique_weighted*100)
        pct_genome = '{:.1f}%'.format(result.f_match*100)

        if is_abundance:
            name = result.match._display_name(screen_width - 41)
            average_abund ='{:.1f}'.format(result.average_abund)
            print_results('{:9}   {:>7} {:>7} {:>9}    {}',
                      format_bp(result.intersect_bp), pct_query, pct_genome,
                      average_abund, name)
        else:
            name = result.match._display_name(screen_width - 31)
            print_results('{:9}   {:>7} {:>7}    {}',
                      format_bp(result.intersect_bp), pct_query, pct_genome,
                      name)

        # write out
        if result_writer is None:
            result_writer = result.init_dictwriter(outfp)
        result.write(result_writer)

        outfp.flush()
        sys.stdout.flush()

    csv_writer.close()
        
    if found:
        # use last result!
        if is_abundance and result:
            p_covered = result.sum_weighted_found / result.total_weighted_hashes
            p_covered *= 100
            print_results(f'the recovered matches hit {p_covered:.1f}% of the abundance-weighted query.')

        print_results(f'the recovered matches hit {sum_f_uniq_found*100:.1f}% of the query k-mers (unweighted).')


if __name__ == '__main__':
    sys.exit(main())
