import operator
from functools import reduce

def delete_keys_from_dict(dictionary, key_list):
    dictionary = reduce(dict.get, key_list[0:-1], dictionary)
    del dictionary[key_list[-1]]

def getFromDict(dataDict, mapList):
    return reduce(operator.getitem, mapList, dataDict)


def setInDict(dataDict, mapList, value):
    """
    sets value to nested dict
    """
    getFromDict(dataDict, mapList[:-1])[mapList[-1]] = value


def appendInDictList(dataDict, mapList, value):
    """
    appends value to nested dict with list-values
    """
    getFromDict(dataDict, mapList[:-1])[mapList[-1]].append(value)


def nested_val_set(dic, keys, value):
    for key in keys[:-1]:
        dic = dic.setdefault(key, {})
    dic[keys[-1]] = value


def nested_dict_pairs_iterator(dict_obj):
    ''' This function accepts a nested dictionary as argument
        and iterate over all values of nested dictionaries
    '''
    # Iterate over all key-value pairs of dict argument
    for key, value in dict_obj.items():
        # Check if value is of dict type
        if isinstance(value, dict):
            # If value is dict then iterate over all its values
            for pair in nested_dict_pairs_iterator(value):
                yield (key, *pair)
        else:
            # If value is not dict type then yield the value
            yield (key, value)


def merge(a, b, path=None):
    "merges b into a"
    if path is None: path = []
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                merge(a[key], b[key], path + [str(key)])
            elif a[key] == b[key]:
                pass # same leaf value
            else:
                raise Exception('Conflict at %s' % '.'.join(path + [str(key)]))
        else:
            a[key] = b[key]
    return a
