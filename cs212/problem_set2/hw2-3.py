# --------------
# User Instructions
#
# Write a function, longest_subpalindrome_slice(text) that takes
# a string as input and returns the i and j indices that
# correspond to the beginning and end indices of the longest
# palindrome in the string.
#
# Grading Notes:
#
# You will only be marked correct if your function runs
# efficiently enough. We will be measuring efficency by counting
# the number of times you access each string. That count must be
# below a certain threshold to be marked correct.
#
# Please do not use regular expressions to solve this quiz!


def is_palindrone(text, length):
    i = 0
    j = length - 1
    while i<j:
        if text[i] != text[j]:
            return False
        i += 1
        j -= 1
    return True

def size(pair):
    return pair[1] - pair[0]

def longest_subpalindrome_slice(text):
    "Return (i, j) such that text[i:j] is the longest palindrome in text."
    if text == '': return (0,0)
    text = text.lower()
    dp = (0,1)
    for i in range(1, len(text)):
        for j in range(0, i+1-size(dp)):
            if is_palindrone(text[j:i+1], i-j+1):
                if i-j+1 > size(dp):
                    dp = (j, i+1)
                    break;

    return dp

def test():
    L = longest_subpalindrome_slice
    assert L('racecar') == (0, 7)
    assert L('Racecar') == (0, 7)
    assert L('RacecarX') == (0, 7)
    assert L('Race carr') == (7, 9)
    assert L('') == (0, 0)
    assert L('something rac e car going') == (8,21)
    assert L('xxxxx') == (0, 5)
    assert L('Mad am I ma dam.') == (0, 15)
    return 'tests pass'


print test()
