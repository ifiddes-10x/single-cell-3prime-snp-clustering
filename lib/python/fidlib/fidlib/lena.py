"""
Convenience wrappers for accessing lena IDs

"""
import json
import urllib2
import glob


def get_lena_path(lena_id):
    """
    Returns the path on marsoc based on a lena ID
    :param lena_id: ID
    :return: string
    """
    url = "http://marsoc/api/shimulate/{}".format(lena_id)

    data = urllib2.urlopen(url).read()

    as_json=json.loads(data)

    fqps = as_json['fastq_paths']

    base = ''
    if len(fqps)  > 1:
        base = "/mnt/analysis/marsoc/pipestances/{}".format(lena_id)
    else:
        x = None
        for y in fqps:
            x = y
        if x is not None:
            base = "/mnt/analysis/marsoc/pipestances/{}".format(x)
    return glob.glob(base + "/*/{}/HEAD".format(lena_id))[0]
