def bin_search(tuples, x):  # dla zmniejszenia asymptotycznej złożoności obliczeniowej
    i, j = 0, len(test) - 1
    while True:
        print(f"i:{i} j:{j}")
        middle_ind = (i+j)//2
        middle_val, middle_x = tuples[middle_ind]
        if middle_x == x:
            return middle_val
        elif i >= j - 1:
            return None
        elif x < middle_x:
            j = middle_ind
        else:
            i = middle_ind

test = [(10,2), (4,5), (7,8), (3, 21), (5, 34)]
print(bin_search(test, 8))