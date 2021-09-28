# RtarD - finding strongly connected components in a digraph using R

This code implements in ***R*** the uber-famous ***Tar***jan's algorithm for finding strongly connected components in a digraph by ***D***epth-first Search

## what's the point?
I hear you. Obviously there are plenty of better-performing, readily-available alternatives out there to work out the strongly connected components (SCC) in a digraph using Depth-First Search (DFS) - an algorithm developed by R. Tarjan in 1972. Examples include built-in functions in Gephi (https://gephi.org/); some Java freebies straight out of the book Algorithms 4th ed (https://algs4.cs.princeton.edu/42digraph/); and the multi-platform package igraph (https://igraph.org/redirect.html). To add insult to injury, the latter is readily available for R.

And yet I am not aware of DFS-SCC functions written for R that are not wrappers for code written in another language. Perhaps it's obvious why not R: the idea of DFS builds on a recursion - a function calling itself. Writing  recursions - as in the original 1972 pseudocode - is not straightforward using the R syntax (unlike say JavaScript). In this code I use "while" loops instead of recursions. That might makes you cringe, but the whole point here is self-study, so my conscience is clean.

Also, in this code the search for SCC is coupled with the problem of finding a block-triangular permutation for the adjacency matrix of the associated digraph (works under certain conditions). I am not aware of this feature being explored jointly with the problem of finding SCC - so I figured that would make some academic geek happy.

## Useful references 
 1) Original academic paper: Page 157, Tarjan (1972) doi:10.1137/0201010
 2) Deo (1974) Graph Theory. ISBN: 9780486807935 [but focuses on undirected graphs]
 3) wikipedia: https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm#The_algorithm_in_pseudocode
 4) Princeton web resources linked to the book "Algorithms 4th ed": https://algs4.cs.princeton.edu/42digraph/
 5) Duff and Reid (1978) doi:10.1145/355780.355785 - FORTRAN code link to block triangularisation. But uses 'go to' statements to jump around...
 6) Strang (1986) Introduction to Applied Mathematics. ISBN 0-9614088-0-4, Ch. 16 (doesn't implement Tarjan but refers to earlier work)
 7) Hume and Plemmons (1981) p 272 in thid Google book https://books.google.co.uk/books?id=pEMsAQAAIAAJ&pg=PA272&lpg=PA272&dq=duff+and+reid+block+triangularization&source=bl&ots=zlTT95Usx-&sig=ACfU3U3GkBU0av0mI459KObNAqyrc-1sTw&hl=en&sa=X&ved=2ahUKEwiWkKnBlKbyAhUSesAKHXmZArcQ6AF6BAgREAM#v=onepage&q=duff%20and%20reid%20block%20triangularization&f=false
