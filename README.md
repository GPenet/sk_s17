# sk_s17
search of 17 clues puzzles in solution grids

This code is a variant of the same action included in the repository skmpp2
This is for the distribution 665 (in bands and stacks). 
Some changes would be needed to search a band 3 having more than 6 clues (blue pass 1)

In the first code, the external double loop combines all valid bands 1 to all valid bands2.
This is derived from blue's findings and is explained in the following thread
http://forum.enjoysudoku.com/scan-solution-grids-for-17-clues-as-of-blue-t34012.html

Here the external loop is in the band with the lowest count (band 1 or band 2).
The second band is then produced in a way similar to the processing of band 3 in the previous code.

The target is to try to speed up the process for bands having a big numbers of valids band1/band2.
First tests lets hope a 40% improvement in run time.

second shot started end September 2019===========================

The first version has been released as v1.0
The second version will be an hybrid process with several key changes
Again, the target is to speed up the process 

The new guidelines are the following

The external loop remains the band (A) (out of bands 1;2) with the lowest count of clues. band B is the other band (1;2)
As before, the first band   is processed by sub lots having the same first 3 clues.
A new set of GUAs is added with UAs in band 3 of size 4 in (2 boxes, 2 rows,4 columns)
Va A permanent set of the smallest 64 UAs is defined for each pair A;B
Vb A permanent set of GUAs sockets is defined with a limit of 256 sockets 
  81 GUAs2 81 GUAs3 plus room for the new sockets of foru columns 
All bands (so including each band 3) has at the start a full 6 clues expansion (3/6 to have the level 3 clues)
Each valid band 3 has its vector Va settled at the start
Each valid band B has its vector Va settled at the beginning of the A;B process

For each 3 first clues A, a reduction of the tables is done (UAs, GUAs ...)
For each 5 clues valid band A
  The count of the minimum number of clues in B is done
    If the count is 6, a direct expansion is done
    if the count is below 6, a vector check is done using the B 6 clues table of valid bands B
    
This gives at the end a set of "valid bands 1+2"
For each valid band 1+2
  Gangster situation is extracted 
  For each band 3, the count of the minimum number of clues is done.
   if the count is 6 a direct expansion is done
   if the count is below 6, a vector check is done using the band 3 specific 6 clues table
 
