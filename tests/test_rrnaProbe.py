# import functions to test
from xpresspipe.rrnaProbe import rrnaProbe, getMatchLeft, getMatchRight, addEntry, combineSeqs

# test getMatchRight (checks for a match on the right side of the new string
# with the left side of any previously added string

#Test 1 (getMatchRight()): positive match on right
#Dummy test data
min_overlap = 5
seq = "supercalif"
freq = 10
combined_seqs = [["califragilisticexpialidocious",12], ["doodle", 2]]

combined_seqs,to_add,found = getMatchRight(seq, freq, combined_seqs, min_overlap)
assert(combined_seqs == [["doodle", 2]])
assert(found)
assert(to_add == ['supercalifragilisticexpialidocious', 22])


#Test 2 (getMatchRight()): negative match
#Dummy test data
min_overlap = 5
seq = "supercali"
freq = 10
combined_seqs = [["califragilisticexpialidocious",12], ["doodle", 2]]

combined_seqs,to_add,found = getMatchRight(seq, freq, combined_seqs, min_overlap)

assert(combined_seqs == [["califragilisticexpialidocious",12], ["doodle", 2]])
assert(not found)
assert(to_add == [])

#Test 3 (getMatchLeft()): positive match on left
#Dummy test data
min_overlap = 5
seq = "expialidocious"
freq = 10
combined_seqs = [["supercalifragilisticexpia",12], ["doodle", 2]]

combined_seqs,to_add,found = getMatchLeft(seq, freq, combined_seqs, min_overlap)

assert(combined_seqs == [["doodle", 2]])
assert(found)
assert(to_add == ["supercalifragilisticexpialidocious",22] )

#Test 4 (getMatchLeft()): negative match
#Dummy test data
min_overlap = 5
seq = "xpialidocious"
freq = 10
combined_seqs = [["supercalifragilisticexpia",12], ["doodle", 2]]

combined_seqs,to_add,found = getMatchLeft(seq, freq, combined_seqs, min_overlap)

assert(combined_seqs == [["supercalifragilisticexpia",12], ["doodle", 2]])
assert(not found)
assert(to_add == [] )

#Test 5 (addEntry()): new seq
#Dummy test data
min_overlap = 5
seq = "supercalifragilisticexpialidocious"
freq = 10
combined_seqs = [["doodle", 2]]

combined_seqs = addEntry([seq,freq], combined_seqs, min_overlap)

assert(combined_seqs == [["doodle", 2], ["supercalifragilisticexpialidocious",10]])

#Test 6 (addEntry()): new seq matching old one
#Dummy test data
min_overlap = 5
seq = "supercalifrag"
freq = 10
combined_seqs = [["doodle", 2], ["califragilisticexpialidocious",10]]

combined_seqs = addEntry([seq,freq], combined_seqs, min_overlap)

assert(combined_seqs == [["doodle", 2], ["supercalifragilisticexpialidocious",20]])

#Test 7 (addEntry()): new seq contained in old
#Dummy test data
min_overlap = 5
seq = "califrag"
freq = 10
combined_seqs = [["doodle", 2], ["supercalifragilisticexpialidocious",10]]

combined_seqs = addEntry([seq,freq], combined_seqs, min_overlap)

assert(combined_seqs == [["doodle", 2], ["supercalifragilisticexpialidocious",20]])

#Test 8 (addEntry()): old seq contained in new
#Dummy test data
min_overlap = 5
seq = "supercalifragilisticexpialidocious"
freq = 10
combined_seqs = [["doodle", 2], ["califrag",10]]

combined_seqs = addEntry([seq,freq], combined_seqs, min_overlap)

assert(combined_seqs == [["doodle", 2], ["supercalifragilisticexpialidocious",20]])

#Test 9 (combineSeqs()):
#Dummy test data
min_overlap = 5
entries = [
    ["doodle", 10],
    ["docious",10],
    ["super",10],
    ["supercalifragilisticexpialidocious", 10]
]

combined_seqs = combineSeqs(entries, min_overlap)

assert(combined_seqs == [["doodle", 10], ["supercalifragilisticexpialidocious", 30]])
