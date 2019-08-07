# sk_s17
search of 17 clues puzzles in solution grids

This code is a variant of the same action included in the repository skmpp2

In the first code, the external double loop combines all valid bands 1 to all valid bands2.
This is derived from blue's findings and is explained in the following thread
http://forum.enjoysudoku.com/scan-solution-grids-for-17-clues-as-of-blue-t34012.html

Here the external loop is in the band with the lowest count (band 1 or band 2).
The second band is then produced in a similar way to the processing of band 3 in the previous coe.

The target is to try to speed up the process for bands having a big numbers of valids band1/band2.
First tests lets hope a 40% improvement in run time.
