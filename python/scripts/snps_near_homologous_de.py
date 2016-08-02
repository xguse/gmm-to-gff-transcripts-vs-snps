"""Get info for those SNPs within a certain distance of genomic features whose homologs were differentially expressed in the external species."""

import click
print = click.echo

import pandas as pd

from multiprocessing import Pool


bed_headers = ["chrom",
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
               "blockStarts"]

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
    # if not isinstance(knearest_paths,list):
    #     knearest_paths =[knearest_paths]

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



def add_url_col(results, url, bed, bed_headers):
    """Add url column to results df."""
    print("begun adding urls.")
    # correct the templating of the url required by snakemake interpreting it as wildcards
    url = url.replace("[chrom]","{chrom}").replace("[left]","{left}").replace("[right]","{right}")

    bed = pd.read_csv(bed, sep='\t',names=bed_headers)[["chrom", "chromStart", "chromEnd", "name"]]
    bed = bed.rename(columns={"name":"gene_id_internal"})

    combo = pd.merge(left=results, right=bed, how='inner',
                     on="gene_id_internal", left_on=None, right_on=None)

    # count_cols = ['chromStart','chromEnd','SNP_start','SNP_end']
    combo['left'] = combo.apply(lambda row: min(row[['chromStart','chromEnd','SNP_start','SNP_end']]), axis=1)
    combo['right'] = combo.apply(lambda row: max(row[['chromStart','chromEnd','SNP_start','SNP_end']]), axis=1)
    combo['url'] = combo.apply(lambda row: url.format(chrom=row["chrom"],left=row["left"],right=row["right"]), axis=1)

    results = pd.merge(left=results, right=combo[["gene_id_internal","url"]], how='inner',
                     on="gene_id_internal", left_on=None, right_on=None)

    return results.drop_duplicates()



def print_arguments(**kwargs):
    """Print a dictionary of the options provided.

    Helpful for making sure the command is being parsed correctly.
    """
    print(kwargs)






@click.command()
@click.option("-d", "--distance", type=int, default=0,
              help="Distance of SNP from feature [default: 0]")
@click.option("-c", "--chunksize", type=int, default=10000,
              help="Rows of diff expr table to read in at one time [default: 10000]")
@click.option("-p", "--parallel", type=click.Choice(['yes', 'no']), default='yes',
              help="Run in parallel [default: yes]")
@click.option("-u", "--url", type=str, default=False,
                help="A 'templated' genome browser url to an ensembl-based genome browser (see above). If provided an extra column will be included in output that includes a link to the chromosome location with range such that the SNP and full coords of the 'gene_id_internal' are included. If this is provided, the --bed option MUST be provided as well.")
@click.option("-b", "--bed", type=click.Path(exists=True),
                help="A path to a bed formated coordinate file for the 'internal' species of interest. (REQUIRED if --url is used).")
@click.option("-n", "--dryrun", is_flag=True,
              help="Do nothing. Just print the options dictionary for sanity")
@click.argument("knearest", type=click.Path(exists=True), nargs=-1, required=True)
@click.argument("diffexpr", type=click.Path(exists=True), nargs=1)
@click.argument("output", type=click.Path(), nargs=1)
def run(distance,chunksize,parallel,url,bed,dryrun,knearest,diffexpr,output):
    """Get info for those SNPs within a certain distance of genomic features whose homologs were differentially expressed in the external species.

    \b
    Args:
        KNEAREST    One or more paths representing a `nearest_features_beds` file (REQIRED).
        DIFFEXPR    A single path representing a `ids_with_diff_expr` file (REQIRED).
        OUTPUT      A single path representing where to write output (REQIRED).

    Write file representing which SNPs are within a certain distance of a feature in `internal`
    species that were diff expressed in `external` species.

    If --url is used the format of the url should include place-holders for the chromosome, and range coordinates similar to the following:

    \b
    https://www.vectorbase.org/Glossina_fuscipes/Location/View?r=[chrom]:[left]-[right]
    """
    # validations
    if url:
        if not bed:
            raise ValueError("You must provide a --bed if you provide a --url.")


    if dryrun:
        # YOU HAVE TO UPDATE THIS MANUALLY WHEN THE ARGS TO MAIN CHANGE!
        print_arguments(distance=distance,
                        chunksize=chunksize,
                        parallel=parallel,
                        dryrun=dryrun,
                        knearest=knearest,
                        diffexpr=diffexpr,
                        output=output,
                        url=url,
                        bed=bed,)
        exit(0)


    print("BEGIN: join_all_k_nearest_with_de")
    results = join_all_k_nearest_with_de(knearest_paths=knearest,
                                         distance=distance,
                                         de_path=diffexpr,
                                         chunksize=chunksize,
                                         parallel=parallel)
    print("END: join_all_k_nearest_with_de")

    if url:
        results = add_url_col(results, url, bed, bed_headers)

    results.to_csv(output, index=False)



if __name__ == '__main__':
    run()
