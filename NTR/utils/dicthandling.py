import operator
from functools import reduce


def getFromDict(dataDict, mapList):
    return reduce(operator.getitem, mapList, dataDict)


def setInDict(dataDict, mapList, value):
    """
    sets value to nested dict
    :param dataDict:
    :param mapList:
    :param value:
    """
    getFromDict(dataDict, mapList[:-1])[mapList[-1]] = value


def appendInDictList(dataDict, mapList, value):
    """
    appends value to nested dict with list-values
    :param dataDict:
    :param mapList:
    :param value:
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


