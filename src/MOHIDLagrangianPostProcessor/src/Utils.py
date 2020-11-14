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

