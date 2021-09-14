This is an R implementation of Tarjan's algorithm for finding Strongly Connected Components (SCC) in a digraph by Depth-First Search (DFS)

Contacts:
   e.settanni@eng.cam.ac.uk
   
Developed from scratch for self-study during summer 2021. it is based on my own reading of the original pseudocode.

useful references 
 1) Original academic paper: Page 157, Tarjan (1972) doi:10.1137/0201010
 2) Deo (1974) Graph Theory [but for undirected graphs]
 3) wikipedia: https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm#The_algorithm_in_pseudocode
 4) Princeton web resources linked to the book "Algorithms 4th ed": https://algs4.cs.princeton.edu/42digraph/
 5) Duff and Reid (1978) doi:10.1145/355780.355785 - FORTRAN code link to block triangularisation. But uses 'go to' statements hence cannot replicate with R
 6) Strang (1986) Introduction to Applied Mathematics, Ch. 16 (doesn't implement Tarjan but refers to earlier work
 7) Hume D., Litsey J., and Plemmons (1981) p 272 - see Google books https://books.google.co.uk/books?id=pEMsAQAAIAAJ&pg=PA272&lpg=PA272&dq=duff+and+reid+block+triangularization&source=bl&ots=zlTT95Usx-&sig=ACfU3U3GkBU0av0mI459KObNAqyrc-1sTw&hl=en&sa=X&ved=2ahUKEwiWkKnBlKbyAhUSesAKHXmZArcQ6AF6BAgREAM#v=onepage&q=duff%20and%20reid%20block%20triangularization&f=false


Benchmarks include BUILT-IN functions in e.g. Gephi - claims it implement Tarjan's algorithm - and igraph.
