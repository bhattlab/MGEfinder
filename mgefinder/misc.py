import warnings
warnings.filterwarnings("ignore")

from scipy.stats import poisson
from bisect import bisect_left, bisect_right


def revcomp(read):
    reversed_seq = ''
    for l in reversed(read.upper()):
        if l == 'A':
            reversed_seq += 'T'
        elif l == 'T':
            reversed_seq += 'A'
        elif l == 'C':
            reversed_seq += 'G'
        elif l == 'G':
            reversed_seq += 'C'
        else:
            reversed_seq += l
    return reversed_seq


def poisson_test_greater(x, mu):

    return 1 - poisson.cdf(x-1, mu)


def takeClosestSmaller(sorted_list, number):

    if len(sorted_list) == 0:
        return None

    pos = bisect_left(sorted_list, number)

    if pos == 0:
        if sorted_list[0] >= number:
            return None
        else:
            return sorted_list[0]

    if pos == len(sorted_list):
        if sorted_list[-1] < number:
            return sorted_list[-1]
        else:
            return sorted_list[-2]

    if sorted_list[pos] < number:
        return sorted_list[pos]
    else:
        return sorted_list[pos-1]


def takeClosestLarger(sorted_list, number):

    if len(sorted_list) == 0:
        return None

    pos = bisect_right(sorted_list, number)

    if pos == len(sorted_list):
        if sorted_list[-1] <= number:
            return None
        else:
            return sorted_list[-1]

    if pos == 0:
        if sorted_list[0] > number:
            return sorted_list[0]
        else:
            return sorted_list[1]

    if sorted_list[pos] > number:
        return sorted_list[pos]
    else:
        return sorted_list[pos+1]



if __name__ == "__main__":
    print(takeClosestSmaller([], 100), 100)
    print()
    print(takeClosestSmaller([1], 0), 0)
    print(takeClosestSmaller([1], 1), 1)
    print(takeClosestSmaller([1], 2), 2)
    print()
    print(takeClosestSmaller([1, 2], 0), 0)
    print(takeClosestSmaller([1, 2], 1), 1)
    print(takeClosestSmaller([1, 2], 2), 2)
    print(takeClosestSmaller([1, 2], 3), 3)
    print()
    print(takeClosestSmaller([1, 2, 3], 0), 0)
    print(takeClosestSmaller([1, 2, 3], 1), 1)
    print(takeClosestSmaller([1, 2, 3], 2), 2)
    print(takeClosestSmaller([1, 2, 3], 3), 3)
    print(takeClosestSmaller([1, 2, 3], 4), 4)
    print(takeClosestSmaller([1, 2, 3], 5), 5)
    print()
    print()
    print(takeClosestLarger([], 100), 100)
    print()
    print(takeClosestLarger([1], 0), 0)
    print(takeClosestLarger([1], 1), 1)
    print(takeClosestLarger([1], 2), 2)
    print()
    print(takeClosestLarger([1, 2], 0), 0)
    print(takeClosestLarger([1, 2], 1), 1)
    print(takeClosestLarger([1, 2], 2), 2)
    print(takeClosestLarger([1, 2], 3), 3)
    print()
    print(takeClosestLarger([1, 2, 3], 0), 0)
    print(takeClosestLarger([1, 2, 3], 1), 1)
    print(takeClosestLarger([1, 2, 3], 2), 2)
    print(takeClosestLarger([1, 2, 3], 3), 3)
    print(takeClosestLarger([1, 2, 3], 4), 4)
    print(takeClosestLarger([1, 2, 3], 5), 5)