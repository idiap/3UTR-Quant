import pandas as pd

def quantile_normalize(df: pd.DataFrame):
    """Performs quantile normalization across the columns.
    This function uses the function normalizeQuantiles from R.

    Args:
        df: pd.DataFrame on which quantile normalization will be performed on the columns.

    Returns:
        pd.DataFrame: The normalized dataframe.
    """

    import copy
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr

    with(robjects.default_converter + pandas2ri.converter).context():
        r_bulk = robjects.conversion.get_conversion().py2rpy(df)

    limma = importr("limma")
    r_bulk_norm = limma.normalizeQuantiles(r_bulk)
    df_norm = pd.DataFrame(
        dict(zip(r_bulk_norm.names, list(r_bulk_norm)))
    )
    df_norm.index = df.index
    return df_norm


def CPM_normalize(df):
    """Performs Count per Million normalization.

    Args:
        df: pd.DataFrame on which CPM normalization will be performed on the columns.

    Returns:
        pd.DataFrame: The normalized dataframe.
    """
    from bioinfokit.analys import norm
    nm = norm()
    nm.cpm(df=df)
    df_norm = nm.cpm_norm
    return df_norm

