"""Get info for those SNPs within a certain distance of genomic features whose homologs were differentially expressed in the external species.

Usage:
  snps_near_homologous_de.py (--knearest=<path>... --diffexpr=<path> --output=<path>) [options]

Options:
  -h --help             Show this screen.
  -v --version          Show version.
  --knearest=<path>     One or more paths representing a `nearest_features_beds` file (REQIRED).
  --diffexpr=<path>     A single path representing a `ids_with_diff_expr` file (REQIRED).
  --output=<path>       A single path representing where to write output (REQIRED).
  -d --distance=<INT>   Distance of SNP from feature [default: 0].
  -c --chunksize=<INT>  Rows of diff expr table to read in at one time [default: 10000].
  -p --parallel         Run in parallel (options are 'yes' or 'no') [default: 'yes'].
  -n --dryrun           Do nothing. Just print the options dictionary for sanity.

Write file representing which SNPs are within a certain distance of a feature in `internal`
species that were diff expressed in `external` species.
"""

from docopt import docopt as opt

import pandas as pd

from multiprocessing import Pool

import ipdb



def load_k_nearest_bed_by_distance(path, distance, keep_cols=None, rename_cols=None):
    """Load and reconfigure the calculated k_nearest file then return it as dataframe."""
    headers = ["SNP_chrom",
           "SNP_start",
           "SNP_end",
           "feature_set_name",
           "chrom",
           "chromStart",
           "chromEnd",
           "name",
           "score",
           "strand",
           "thickStart",
           "thickEnd",
           "itemRgb",
           "blockCount",
           "blockSizes",
           "blockStarts",
           "distance"
          ]
    k_nearest = pd.read_csv(path, sep="\t", names=headers)

    filtered_by_d = k_nearest.query(""" abs(distance) <= {distance} """.format(distance=distance))

    if rename_cols is not None:
        assert isinstance(rename_cols,dict)
        filtered_by_d = filtered_by_d.rename(columns=rename_cols).copy()

    if keep_cols is not None:
        assert isinstance(keep_cols,list)
        return filtered_by_d[keep_cols]
    else:
        return filtered_by_d

def load_de_genes_tx(path, chunksize=1000):
    """Load the compiled differential expression file by chunks and return iterator of dataframes."""
    de = pd.read_csv(path, chunksize=chunksize)
    return de


def join_all_k_nearest_with_de(knearest_paths, distance, de_path, chunksize, parallel):
    """Return dataframe where SNPs are within a certain distance of a feature in `internal` species that were diff expressed in `external` species."""
    if not isinstance(knearest_paths,list):
        knearest_paths =[knearest_paths]

    if parallel == 'yes':
        job_args = []

        # Generate each job's arguments
        for knearest_path in knearest_paths:
            job_args.append([knearest_path, distance, de_path, chunksize])

        # Run the jobs in parallel with the workers pool
        with Pool(8) as workers:
            results_tmp = workers.starmap(join_one_k_nearest_with_de, job_args)

        results = []
        for r in results_tmp:
            results.extend(r)
    elif parallel == 'no':
        results =[]
        for path in knearest_paths:
            results_tmp = join_one_k_nearest_with_de(knearest_path=path, distance=distance, de_path=de_path, chunksize=chunksize)
            results.extend(results_tmp)
    else:
        raise ValueError("Valid values for '--parallel' are ['yes','no']. You provided: {value}".format(value=parallel))

    return pd.concat(results)


def join_one_k_nearest_with_de(knearest_path, distance, de_path, chunksize):
    """Perform merge of `proximal_id`s for each `feature_set_name` for a single k_nearest file and return a list of the merged chunks."""
    #

    print("joining {knearest} against {de_path}\n".format(knearest=knearest_path.split('/')[-1],de_path=de_path.split('/')[-1]))


    kn_df = load_k_nearest_bed_by_distance(path=knearest_path,
                                           distance=distance,
                                           keep_cols=["SNP_chrom","SNP_start","SNP_end","feature_set_name","proximal_id"],
                                           rename_cols={"name": "proximal_id"})

    kn_by_feature_set = kn_df.groupby('feature_set_name')

    results = []
    for chunk in load_de_genes_tx(path=de_path, chunksize=chunksize):

        for feature_set_name, group in kn_by_feature_set:

            if feature_set_name == 'official_annotations':
                join_column = 'gene_id_internal'
            elif feature_set_name == 'novel_mapped_tx':
                join_column = 'tcons_id'

            merged_chunk = group.merge(right=chunk, how='inner',
                                       on=None, left_on="proximal_id", right_on=join_column,
                                       left_index=False, right_index=False)

            results.append(merged_chunk.copy())

    # results_flattened = [result for sublist in results for result in sublist]
    return results




def main(nearest_features_beds, snp_distance_from_gene, ids_with_diff_expr, de_table_chunksize, run_parallel):
    """Get info for those SNPs within a certain distance of genomic features whose homologs were differentially expressed in the external species.

    Write file representing which SNPs are within a certain distance of a feature in `internal` species that were diff expressed in `external` species.
    """
    results = join_all_k_nearest_with_de(knearest_paths=nearest_features_beds,
                                         distance=snp_distance_from_gene,
                                         de_path=ids_with_diff_expr,
                                         chunksize=de_table_chunksize,
                                         parallel=run_parallel)

    results.to_csv(snps_near_homologous_de_path, index=False)



if __name__ == '__main__':
    args = opt(__doc__, version='2016.04.05.1')

    # ipdb.set_trace()
    if args["--dryrun"]:
        print(args)
        exit(0)


    # Settings
    snp_distance_from_gene = args["--distance"]
    assert isinstance(snp_distance_from_gene, int)

    de_table_chunksize = args["--chunksize"]
    assert isinstance(de_table_chunksize, int)

    run_parallel = args["--parallel"]


    # input
    nearest_features_beds = args["--knearest"]
    ids_with_diff_expr = args["--diffexpr"]


    #output
    snps_near_homologous_de_path = args["--output"]


    main(nearest_features_beds,
         snp_distance_from_gene,
         ids_with_diff_expr,
         de_table_chunksize,
         run_parallel)
