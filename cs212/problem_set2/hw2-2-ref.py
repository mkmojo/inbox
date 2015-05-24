import itertools

def floor_puzzle():
    floors = bottom, _, _, _, top = [1, 2, 3, 4, 5]
    orderings = list(itertools.permutations(floors))
    for (Hopper, Kay, Liskov, Perlis, Ritchie) in orderings:
        if (Kay != bottom
            and Liskov != top  and Liskov != bottom
            and Perlis > Kay
            and abs(Ritchie - Liskov) > 1
            and abs(Kay - Liskov) >1 ):
            return [Hopper, Kay, Liskov, Perlis, Ritchie]
print floor_puzzle()
