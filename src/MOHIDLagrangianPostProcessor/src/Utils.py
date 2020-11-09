# -*- coding: utf-8 -*-


def filterByTrue(_list, mask):
    return [_list[i] for i in range(0, len(_list)) if mask[i]]


def isWordinListElements(word, _list):
    return [word in _listElement for _listElement in _list]


def groupByWords(_list, words):
    for word in words:
        mask = isWordinListElements(word, _list)
        _list = filterByTrue(_list, mask)
    return _list


def getTotalAmmountEmitted(dataset):
    if 'n_counts_global' in dataset:
        n = dataset['n_counts_global'].diff('time').sum().values
        return n


def getMeanStdGlobal(dataset, dim=None):
    if dim is None:
        return dataset['n_counts_global'].mean(), dataset['n_counts_global'].std()
    else:
        return dataset['n_counts_global'].mean(dim=dim), dataset['n_counts_global'].std(dim=dim)


def getMaxStdGlobal(dataset):
    return dataset['n_counts_global'].max()

