## expand function from snakemake
import itertools
import string
import collections

def expand(*args, **wildcards):
    """
    This function was adopted from snakemake
    Works exactly like snakemake's expand function

    Expand wildcards in given filepatterns.
    Arguments
    *args -- first arg: filepatterns as list or one single filepattern,
        second arg (optional): a function to combine wildcard values
        (itertools.product per default)
    **wildcards -- the wildcards as keyword arguments
        with their values as lists. If allow_missing=True is included
        wildcards in filepattern without values will stay unformatted.
    """
    filepatterns = args[0]
    if len(args) == 1:
        combinator = itertools.product
    elif len(args) == 2:
        combinator = args[1]
    if isinstance(filepatterns, str) or isinstance(filepatterns, Path):
        filepatterns = [filepatterns]

    def path_to_str(f):
        if isinstance(f, Path):
            return str(f)
        return f

    filepatterns = list(map(path_to_str, filepatterns))

    if any(map(lambda f: getattr(f, "flags", {}), filepatterns)):
        raise WorkflowError(
            "Flags in file patterns given to expand() are invalid. "
            "Flags (e.g. temp(), directory()) have to be applied outside "
            "of expand (e.g. 'temp(expand(\"plots/{sample}.pdf\", sample=SAMPLES))')."
        )

    # check if remove missing is provided
    format_dict = dict
    if "allow_missing" in wildcards and wildcards["allow_missing"] is True:

        class FormatDict(dict):
            def __missing__(self, key):
                return "{" + key + "}"

        format_dict = FormatDict
        # check that remove missing is not a wildcard in the filepatterns
        for filepattern in filepatterns:
            if "allow_missing" in re.findall(r"{([^}\.[!:]+)", filepattern):
                format_dict = dict
                break

    # remove unused wildcards to avoid duplicate filepatterns
    wildcards = {
        filepattern: {
            k: v
            for k, v in wildcards.items()
            if k in re.findall(r"{([^}\.[!:]+)", filepattern)
        }
        for filepattern in filepatterns
    }

    def flatten(wildcards):
        for wildcard, values in wildcards.items():
            if isinstance(values, str) or not isinstance(
                values, collections.abc.Iterable
            ):
                values = [values]
            yield [(wildcard, value) for value in values]

    formatter = string.Formatter()
    try:
        return [
            formatter.vformat(filepattern, (), comb)
            for filepattern in filepatterns
            for comb in map(format_dict, combinator(*flatten(wildcards[filepattern])))
        ]
    except KeyError as e:
        raise WildcardError("No values given for wildcard {}.".format(e))