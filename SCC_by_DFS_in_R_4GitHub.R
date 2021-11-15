###  Find Strongly Connected Components in a digraph by Depth-First Search [Tarjan]

#### -1.0 - Contacts ####
#
# Author:
#   Ettore Settanni
#   e.settanni@eng.cam.ac.uk
#   Cambridge
#

#### -1.1 - Revision notes ####
# Coding timeline
#   Start:       13 07 2021
#   End:         06 08 2021
#   Fixes:       19 08 2021
#   for GitHub:  14 09 2021
#   REV01        02 11 2021    # loops the lowlink update in backtracking until the node's lowlink value pick up the smallest lowlink amongst its successors
#   REV02        10 11 2021    # a node's successor may haves successors but these may be sinks - in which case, don't update lowlink
#   REV03 (this) 13 11 2021    # fixes how the sack is popped (REV02 was a patch, but would fail). Keeps the intuition of REV01: the iteration is still needed

#   One thing I didn't notice when using 'while loops' is the need for extra iterations when we're done backtracking to ensure a correct lowlink update (affecting stack-popping) 
#   Iterations may be needed before all nodes' lowlink values are actually set to correspond to the lowest lowlink value among their successors.
#   Also until REV02 I've got the stack-popping routine misplaced, which would fail in certain instances.


#### -1.2 - motivation, references & caveats ####

# i have developed this implementation from scratch for self-learning/self-study
# it is based on my own (limited) understadinding of the original pseudocode, and build on flowcharts I sketched during Summer 2021

## PSEUDOCODES
# 1) Original academic paper: Page 157, Tarjan (1972) doi:10.1137/0201010
# 2) Deo (1974) Graph Theory [but for undirected graphs]
# 3) wikipedia: https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm#The_algorithm_in_pseudocode
# 4) Princeton web resources linked to the book "Algorithms 4th ed": https://algs4.cs.princeton.edu/42digraph/

## Link with Block triangular permutation
# 5) Duff and Reid (1978) - FORTRAN code link to block triangularisation. But uses 'go to' statements hence cannot replicate with R
# 6) Strang G. (1986) Introduction to Applied Mathematics, Ch. 16 (doesn't implement Tarjan but refers to earlier work)

## Benchmarks - existing BUILT-IN functions 
# Gephi - claims it implement Tarjan's algorithm
# igraph - obviously has a fucntion for SCC, not sure what runds underneath (C wrapper)


#### 00.1 - initialise: Clear ####
rm(list = ls())
gc()




#### 00.2 - initialise: data ####
# 
# test_m <- matrix(c(
#   0, 1, 0, 0, 0, 0, 0, 0,
#   0, 0, 1, 0, 0, 0, 0, 0,
#   1, 0, 0, 0, 0, 0, 0, 0,
#   0, 1, 1, 0, 1, 0, 0, 0,
#   0, 0, 0, 1, 0, 1, 0, 0,
#   0, 0, 1, 0, 0, 0, 1, 0,
#   0, 0, 0, 0, 0, 1, 0, 0,
#   0, 0, 0, 0, 1, 0, 1, 1
# ),ncol = 8, byrow = TRUE)

## STRANG EXAMPLE intro to applied math p 637
# test_m <- matrix(c(
#   0, 0, 1, 0,
#   1, 0, 0, 1,
#   1, 0, 0, 0,
#   0, 0, 1, 0
# ), ncol = 4, byrow = TRUE)

## Princeton example https://algs4.cs.princeton.edu/42digraph/
# test_m <- matrix(c(
#   0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
#   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#   1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#   0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
#   0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#   0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
#   1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,
#   0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0,
#   0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0,
#   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0,
#   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
#   0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
#   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0
# ),ncol = 13, byrow = TRUE)


#### 00.3 - initialise: automated node labelling ####
n_nodes <- ncol(test_m)
node_names <- formatC(0:(ncol(test_m)-1),width = 2, flag = "0")          # thread: https://stackoverflow.com/questions/8266915/format-number-as-fixed-width-with-leading-zeros
colnames(test_m) <- node_names
rownames(test_m) <- node_names


#### 01.1 - Lookups: create list of successors for each node ####

successors <- lapply(node_names, function(x){
  current_node <- x
  # search the row corresponding to the current node for successors
  all_successors_idx <- expand.grid(current_node, node_names, width = 2, flag = "0")                 # expand on current node to find successors
  all_successors_links <- apply(all_successors_idx, 1, function(y){
    test_m[y[1], y[2]]
  })
  if(length(as.character(all_successors_idx[which(all_successors_links != 0),2])) > 0){
    as.character(all_successors_idx[which(all_successors_links != 0),2])
  } else {
    NA
  }
})
names(successors) <- node_names 

#### 01.1 - Lookups: utils table (for depth descent and backtracking) ####
util_table <- lapply(1:n_nodes, function(z){
  x <- successors[z]
  a <- match(names(x), node_names)
  b <- match(unlist(x), node_names)
  temp_col <- expand.grid(a,b)
  tab_head <- c("i", 
                "j", 
                "level", 
                "v_node", 
                #"v_label", 
                "v_number", 
                "v_lowlink", 
                "successors", 
                "k_successor_idx", 
                "successor_k_is_sink", # (REV02)
                "w_k_node", 
                #"w_k_label", 
                "w_k_number", 
                "w_k_lowlink")
  temp_mat <- matrix(NA, ncol = length(tab_head), nrow = nrow(temp_col))
  colnames(temp_mat) <- tab_head
  rownames(temp_mat) <- rep(names(x),nrow(temp_mat))
  temp_mat[,"v_node"] <- temp_col[,1]
  temp_mat[, "w_k_node"] <- temp_col[,2]
  
  return(temp_mat)
})
util_table_df <- do.call(rbind, util_table)




#### 01.2 - initialise stuff ####

# initialise counter
n <- n_nodes
i <- 0
comp_count <- 0                                                                                                     # component counter

# initialise condition flags
nodes_not_numbered_yet <- n
test0 <- TRUE




# initialise lists
Stack_S <- list()
Component_list <- list()

# initialise table of node numbers
node_numbering <- cbind(rep(NA,length(node_names)),rep(NA,length(node_names)),rep(NA,length(node_names)))
colnames(node_numbering) <- c("node_number","node_lowlink","node_onStack")
rownames(node_numbering) <- node_names


#### 02.1 - Main recursion  ####

# Loop #0: restart at each strongly connected component
dummy4 <- TRUE
while(dummy4){
  level <- 1                                                                                               # I AM NOT ENTIRELY SURE ABOUT THIS SYSTEM I CAME UP WITH. it's like a thread that pulls you back up once the depth-first search reached an end. 
  j <- 0
  # pick the first node not numbered yet
  nodes_to_explore <- which(is.na(node_numbering[, "node_number"]))
  if(nodes_not_numbered_yet != 0){                                                                         # there are nodes not numbered yet
    v <- as.numeric(nodes_to_explore[1])                                                                   # just picking the first unnumbered node. In the first iteration this will be node 1
    v_label <- node_names[v]
  }
  
  # REV01
  flag_extra_round <- TRUE                                                                                 # for use later during backtracking to add an extra round and avoid some issues with updating values at the very end
  extra_round_counter <- 0
  
  # Loop #1: explores successors "depth first", jumping between nodes
  dummy0 <- TRUE
  while(dummy0){
    i <- i + 1                                                                                             # labels nodes as we visit them
    j <- j + 1                                                                                             # position in stack of nodes
    # update node numbering
    Stack_S[j] <- v
    node_numbering[v,"node_number"] <- i
    node_numbering[v,"node_lowlink"] <- i
    node_numbering[v,"node_onStack"] <- 1
    # successors
    w_labels <- unlist(successors[v])
    
    # update records were v features as successor
    back_idx <- which(util_table_df[,"w_k_node"] == v)                                                     
    if (length(back_idx) > 0){  
      util_table_df[back_idx, "w_k_number"] <- node_numbering[v,"node_number"]
      util_table_df[back_idx, "w_k_lowlink"] <- node_numbering[v,"node_lowlink"]
    }
    # check if successor is a sink node
    sink_test <- which(is.na(w_labels))                                                                    
    if (length(sink_test) == 0){
      w <- match(w_labels, node_names)
      n_successors <- length(w)
    } else {
      w <- NA
      n_successors <- 0
    }
    # initialise a progressive counter for incident nodes reachable from the current node
    k <- 0
    
    # Loop #2: explores incident nodes "sequentially" for a given node
    dummy1 <- TRUE
    while(dummy1){
      if(length(sink_test) == 0){                                                                          # the current node has a successor
        k <- k + 1
        w_k <- w[k]                                                                                        # next incident node (neighbor)
        # update main tableau (df version)
        util_table_df_subset <- which(rownames(util_table_df) == v_label)                                  # filter for the current node v
        if(length(util_table_df_subset) > 1){
          if(k == 1){
            util_table_df[util_table_df_subset,][k,"i"] <- i
            util_table_df[util_table_df_subset,][k,"j"] <- j
            util_table_df[util_table_df_subset,][k,"level"] <- level
            util_table_df[util_table_df_subset,][k,"v_node"] <- v
            util_table_df[util_table_df_subset,][k,"v_number"] <- node_numbering[v,"node_number"]
            util_table_df[util_table_df_subset,][k,"v_lowlink"] <- node_numbering[v,"node_lowlink"]
            util_table_df[util_table_df_subset,][k,"successors"] <- n_successors
          }
          util_table_df[util_table_df_subset,][k,"k_successor_idx"] <- k
          # util_table_df[util_table_df_subset,][k,"w_k_node"] <- w_k                                        # could be removed... 
          util_table_df[util_table_df_subset,][k,"w_k_number"] <- node_numbering[w_k,"node_number"]
          util_table_df[util_table_df_subset,][k,"w_k_lowlink"] <- node_numbering[w_k, "node_lowlink"]
        } else {
          if(k == 1){
            util_table_df[util_table_df_subset,"i"] <- i
            util_table_df[util_table_df_subset,"j"] <- j
            util_table_df[util_table_df_subset,"level"] <- level
            util_table_df[util_table_df_subset,"v_node"] <- v
            util_table_df[util_table_df_subset,"v_number"] <- node_numbering[v,"node_number"]
            util_table_df[util_table_df_subset,"v_lowlink"] <- node_numbering[v,"node_lowlink"]
            util_table_df[util_table_df_subset,"successors"] <- n_successors
          }
          util_table_df[util_table_df_subset,"k_successor_idx"] <- k
          # util_table_df[util_table_df_subset,"w_k_node"] <- w_k                                            # could be removed... 
          util_table_df[util_table_df_subset,"w_k_number"] <- node_numbering[w_k,"node_number"]
          util_table_df[util_table_df_subset,"w_k_lowlink"] <- node_numbering[w_k, "node_lowlink"]
        }
        # tests
        test0 <- is.na(node_numbering[w_k,"node_number"])                                                  # we can jump onto this node (depth first)
        test1 <- node_numbering[w_k,"node_number"] < node_numbering[v,"node_number"]                       # the next node has been visited already?
        test2 <- node_numbering[w_k,"node_onStack"] == 1                                                   # the next node is in the stack already OR has been on the stack
      } else {                                                                                             # there is no successor to the current node (sink)
        # tests
        test0 <- FALSE                                                                                     # forces to backtrack if there is no successor
        test1 <- FALSE
        test2 <- FALSE
        # update main tableau (df version)
        util_table_df_subset <- which(rownames(util_table_df) == v_label)                                  # filter for the current node v
        util_table_df[util_table_df_subset,"i"] <- i
        util_table_df[util_table_df_subset,"j"] <- j
        util_table_df[util_table_df_subset,"level"] <- level
        util_table_df[util_table_df_subset,"v_node"] <- v
        util_table_df[util_table_df_subset,"v_number"] <- node_numbering[v,"node_number"]
        util_table_df[util_table_df_subset,"v_lowlink"] <- node_numbering[v,"node_lowlink"]
        util_table_df[util_table_df_subset,"successors"] <- n_successors
        util_table_df[util_table_df_subset,"k_successor_idx"] <- k
        util_table_df[util_table_df_subset,"w_k_node"] <- (-1)
        util_table_df[util_table_df_subset,"w_k_number"] <- (-1)
        util_table_df[util_table_df_subset,"w_k_lowlink"] <- (-1)
      }
      
      test_POP <- FALSE
  
      # THIS IS WHERE WE DECIDE WHETHER TO CONTINUE DEPTH-FIRST, LOOK INTO A NEIGHBOUR, OR BACKTRACK
      if (test0){
        # continue depth-first search
        dummy1 <- FALSE                                                                                    # break Loop 2 (exit sequential exploration of successors for a given node)
        v <- w_k                                                                                           # move on to successor
        v_label <- node_names[v]
        # the below avoids having more than one line with the same level
        level_max <- max(util_table_df[,"level"], na.rm = T)                    
        if (level < level_max){                                                                            # probably we backtracked and re-winded the level counter.This should get us back where we were before the rewind 
          level <- level_max + 1
        } else {
          level <- level + 1
        }
      } else {
        # start "back tracking"
        
        # Loop #3: are there neighbors left to explore?
        dummy5 <- TRUE
        while(dummy5){
          if(k != 0 & k < n_successors) {                                                                  # there are other adjacent nodes left to explore
            dummy5 <- FALSE
            
            # update lowlink values - FIRST TYPE (going downwards)
            if (test1 == TRUE & test2 == TRUE){                                                             # the current k successor has been visited PRIOR to v (has lower number), and it is on the stack already
              node_numbering[v,"node_lowlink"] <- min(node_numbering[v,"node_lowlink"],node_numbering[w_k,"node_number"])
              util_table_df_subset <- which(rownames(util_table_df) == v_label)                              # filter for the current node v
              if(length(util_table_df_subset) > 1){
                util_table_df[util_table_df_subset,][1,"v_lowlink"] <- node_numbering[v,"node_lowlink"]      # always update the first line, corresponding to v
              } else {
                util_table_df[util_table_df_subset,"v_lowlink"] <- node_numbering[v,"node_lowlink"]
              }
              # update lowlink elsewhere
              other_idx <- which(util_table_df[,"w_k_node"] == v)
              if(length(other_idx)>0){
                util_table_df[other_idx, "w_k_lowlink"] <- node_numbering[v,"node_lowlink"]
              }
            }
            
            
          } else {
            dummy3 <- TRUE 
            # Loop #4: revisit all the neighbours already visited
            while (dummy3) {
              if(n_successors !=0 & k !=0){                         # we are not at a sink node / we have not re-visited all the neighbors 
                
                # re-do test
                test0 <- is.na(node_numbering[w_k,"node_number"])                                 # we can jump onto this node (depth first)
                test1 <- node_numbering[w_k,"node_number"] < node_numbering[v,"node_number"]      # the next node has been visited already?
                test2 <- node_numbering[w_k,"node_onStack"] == 1                                  # the next node is in the stack already
                
                
                if (test2){                                                                       # the next node is STILL on the stack PLEASE DO NOT ADD TEST 1 OR WILL MESS UP
                  # update lowlink values - SECOND TYPE (ascending)
                  old_v_lowlink <- node_numbering[v,"node_lowlink"]
                  node_numbering[v,"node_lowlink"] <- min(node_numbering[v,"node_lowlink"],node_numbering[w_k,"node_lowlink"])
                  util_table_df_subset <- which(rownames(util_table_df) == v_label)                              # filter for the current node v
                  if(length(util_table_df_subset) > 1){
                    util_table_df[util_table_df_subset,][1,"v_lowlink"] <- node_numbering[v,"node_lowlink"]      # always update the first line, corresponding to v
                  } else {
                    util_table_df[util_table_df_subset,"v_lowlink"] <- node_numbering[v,"node_lowlink"]
                  }
                  # update lowlink elsewhere
                  other_idx <- which(util_table_df[,"w_k_node"] == v)                                            
                  if(length(other_idx)>0){
                    util_table_df[other_idx, "w_k_lowlink"] <- node_numbering[v,"node_lowlink"]
                  }
                }
                
                
                # go back to the previous neighbour already visited
                k <- k - 1                          
                if (k !=0){
                  w_k <- w[k] 
                  #re-do the test
                  test0 <- is.na(node_numbering[w_k,"node_number"])                                 # we can jump onto this node (depth first)
                  test1 <- node_numbering[w_k,"node_number"] < node_numbering[v,"node_number"]      # the next node has been visited already?
                  test2 <- node_numbering[w_k,"node_onStack"] == 1                                  # the next node is in the stack already
                } else {
                  test2 <- FALSE
                  test1 <- FALSE
                  test0 <- FALSE
                  w_k <- (-1)
                }
              } else {
                # either we are at a sink node or we have re-visited all the neighbors
                
                # NEW: WHEN WE ARE HERE, WE ARE FORCING THE BACKTRACKING - THERE WAS NO PATH LEADING US back to previously encountered nodes
                # therefore we MUST NOT UPDATE LOWLINK when we encounter a node we visited before
                # if v is a root node (lowlink = number), we can pop the stack at this stage
                test_POP <- node_numbering[v,"node_number"] == node_numbering[v,"node_lowlink"]       # from Tarjan's original paper
                
                # BEFORE POPPING
                abort_POP <- (test_POP & level == 1)
                if (abort_POP & flag_extra_round){
                  # we reached the last level (at least until we restart with level 1); WE NEED TO make sure that all nodes' lowlink values matches the lowest  lowlink of their successors.
                  # It happens that the root node gets popped out then the lowlink update occurs, and a bunch of nodes are left without root and can't pop next
                  # if that's not the case it may be that we need more rounds of update
                  
                  extra_round_counter <- extra_round_counter + 1
                  # the below checks that we DON'T HAVE NODES WITH THE LATEST SUCCESSOR'S LOWLINK VALUE BEING lower THAN THE NODE'S LOWLINK
                  # if a node's LOWLINK is larger than the lowest among its successors' LOWLINK the we need to keep updating
                  update_needed <- lapply(Stack_S, function(x){
                    v_lowlink <- node_numbering[x,"node_lowlink"]
                    v_neigb_lowest_lowlink <- min(util_table_df[which(util_table_df[,"v_node"] == x) , "w_k_lowlink"])
                    if(v_neigb_lowest_lowlink > 0){                                                       # sink nodes by convention get "-1" in this field which might cause looping forever
                      if(v_lowlink != v_neigb_lowest_lowlink){
                        1
                      } else {0}
                    } else {0}
                    
                  })
                  updates_test <- sum(unlist(update_needed))
                  if(updates_test > 0 & extra_round_counter < 4){
                    # reset the level to a max but add one and let the loop later pick which level is righ
                    level <- max(util_table_df[which(!is.na(util_table_df[,"level"])),"level"]) + 1          # start all over again, hoping we get all the lowlinks updated
                    test_POP <- FALSE                                                                        # ABORT POPPING FOR NOW
                  } else {
                    flag_extra_round <- FALSE                                                             # break out of this
                    abort_POP <- FALSE
                  }
                  
                  
                }
                
                if (test_POP){
                  # THIS IS WHERE WE POP THE STACK AND POPULATE THE COMPONENT
                  # update count of nodes not numbered yet
                  nodes_not_numbered_yet <- length(which(is.na(node_numbering[,"node_number"])))
                  
                  # NEW: we want to pop the stack but only using v as a root (leave other nodes alone, if any) - different test to enter the stack-popping stage
                  temp_idx <- as.numeric(which(!is.na(node_numbering[,"node_number"])))
                  temp_subset_a <- which(node_numbering[temp_idx, "node_number"] >= node_numbering[v,"node_number"])
                  nodes_visited_so_far <- match(names(temp_subset_a[as.numeric(which(node_numbering[names(temp_subset_a), "node_onStack"] == 1))]), node_names)
                  
                  Stack_S_BACKUP <- Stack_S
                  Component_list_BACKUP <- Component_list
                  # Start of "Pop nodes" sub, but only if we're done with depth-first search
                  
                  # NEW: notice that the loop below should go from node "v" onward
                  for(r in nodes_visited_so_far){
                    x <- as.numeric(node_numbering[r,])
                    if (x[1] == x[2]){                         # the node is a root node
                      comp_count <- comp_count + 1                               # component number
                      # put nodes in component
                      temp_component_list <- lapply(1:length(Stack_S), function(s){
                        y <- Stack_S[[s]]
                        pop_test1 <- node_numbering[y, "node_number"] >= x[1]
                        pop_test2 <- node_numbering[y, "node_lowlink"] == x[2]
                        if(pop_test1 & pop_test2){
                          y
                        }
                      })
                      to_delete <- unlist(temp_component_list)
                      Component_list[[comp_count]] <- to_delete
                      # pop elements out of stack
                      #idx_delete <- which(Stack_S == to_delete)
                      temp_idx_delete <- rep(seq_along(Stack_S), sapply(Stack_S, length))         # thread: https://stackoverflow.com/questions/11002391/fast-way-of-getting-index-of-match-in-list
                      idx_delete <- temp_idx_delete[match(to_delete, unlist(Stack_S))]            
                      Stack_S[idx_delete] <- NULL
                      
                      # update j or the size of the stack will change later
                      j <- j - length(idx_delete)
                      
                      # update table
                      node_numbering[to_delete,"node_onStack"] <- (-1)    # mark elements that have been on stack, but no longer are
                    }
                  }
                }
                
                
                
                dummy3 <- FALSE                                                                                                    # break loop 4 after moving one level up - REGARDLESS of popping the stack
                
                # NEW: modified version of finding the "right level" to go back to
                
                # ignore if we're forcing one last iteration to update lowling - special case
                test_skip <- FALSE
                while(!test_skip & level > 0){
                  level <- level - 1                                                                                             # if there is no successor OR the node at this level is not on the stack, keep going backwards
                  if(level != 0){
                    v_temp <- util_table_df[which(util_table_df[,"level"] == level),"v_node"]
                    level_on_stack <- (node_numbering[v_temp,"node_onStack"] == 1)                                               # skip levels that correspond to nodes that have been removed from stack
                    test_any_successor <- length(which(util_table_df[,"level"] == level & util_table_df[, "successors"] !=0))    # (no successor going up one level)
                    test_skip <- (level_on_stack & test_any_successor > 0)
                  }
                }
                
                # GET BACK TO WHERE YOU LEFT THINGS AT THE PREVIOUS LEVEL
                if (level > 0){
                  back_idx <- which(util_table_df[,"level"] == level & util_table_df[, "successors"] !=0)       # there may be more than one node with the same level if E.G. ONE NEIGHBOUR IS A SINK but the other isn't. We ignore the sink
                  v <- util_table_df[back_idx, "v_node"]                                                  
                  v_label <-  node_names[v]
                  n_successors <- util_table_df[back_idx, "successors"]
                  k_back_idx_a <- as.numeric(which(util_table_df[, "v_node"] == v))
                  util_table_df_subset <- util_table_df[k_back_idx_a,]
                  if(length(k_back_idx_a) > 1){
                    k_back_idx_b <- max(which(!is.na(util_table_df_subset[, "k_successor_idx"])))               # retrieve last successor index explored for the current node
                    k <- util_table_df_subset[k_back_idx_b , "k_successor_idx"]
                    w_k <- util_table_df_subset[k_back_idx_b , "w_k_node"]
                  } else {
                    k <- as.numeric(util_table_df_subset["k_successor_idx"])
                    w_k <- as.numeric(util_table_df_subset["w_k_node"])
                  }
                  w <- as.numeric(util_table_df[k_back_idx_a, "w_k_node"])
                  w_labels <- node_names[w]
                  sink_test <- which(is.na(w_labels))                                                         # check if the successor is a sink node
                  # re-do test
                  test0 <- is.na(node_numbering[w_k,"node_number"])                                 # we can jump onto this node (depth first)
                  test1 <- node_numbering[w_k,"node_number"] < node_numbering[v,"node_number"]      # the next node has been visited already?
                  test2 <- node_numbering[w_k,"node_onStack"] == 1                                  # the next node is in the stack already
                }
              }
            } # end of Loop 4 (dummy 3)
            if(level == 0) {
              dummy5 <- FALSE                                                                                # break loop 3
              dummy1 <- FALSE                                                                                # break loop 2 (ends up in the same place as if there were no successors)
            }
          }  
        }  # end of Loop 3 (dummy 5)
      } #  end of ELSE - "back tracking" sub
    } # end of Loop 2 (dummy1)
    nodes_not_numbered_yet <- length(which(is.na(node_numbering[,"node_number"])))
    nodes_visited_so_far <- as.numeric(which(!is.na(node_numbering[,"node_number"]) & node_numbering[,"node_onStack"] == 1))
    if(!test0 & test_POP){
      if(nodes_not_numbered_yet == 0){                                                                     # break Loop 1 (depths first exploration of successors, jumping between nodes) and move up one level
        dummy4 <- FALSE                                                                                    # break Loop 5: YOU'RE DONE
      }
      dummy0 <- FALSE                                                                                      # break Loop 1
      # if we are here, we'll restart the level counter. re-label the LEVELS used so that they do not interefere with the next component
      util_table_df_BACKUP <- util_table_df
      idx_change_level <- which(!is.na(util_table_df[,"level"]) & util_table_df[,"level"] > 0)
      util_table_df[idx_change_level,"level"] <- (-1)*util_table_df[idx_change_level,"level"]
    }
  } # end of Loop 1 (dummy0)
} # end of Loop 0 (dummy4)

#### 02.2 - Block triangular form (bonus) ####

# Based on Plain English description from: Hume D., Litsey J., and Plemmons (1981) p 272 - see this google book https://books.google.co.uk/books?id=pEMsAQAAIAAJ&pg=PA272&lpg=PA272&dq=duff+and+reid+block+triangularization&source=bl&ots=zlTT95Usx-&sig=ACfU3U3GkBU0av0mI459KObNAqyrc-1sTw&hl=en&sa=X&ved=2ahUKEwiWkKnBlKbyAhUSesAKHXmZArcQ6AF6BAgREAM#v=onepage&q=duff%20and%20reid%20block%20triangularization&f=false
# - all vertices within a strongly connected component are numbered sequentially
# - if there exists and edge from a vertex in component i to a vertex in component j, then all vertices in component i are labelled before all those in component j
n_components <- length(Component_list)
x_components_idx <- expand.grid(1:n_components, 1:n_components)                                                     # equivalent to nested for
x_components_idx <- x_components_idx[, ncol(x_components_idx):1]                                                    # flip
x_comp_link <- apply(x_components_idx, 1, function(x){
  if(x[1] != x[2]){
    comp_nodes_a <- Component_list[[as.numeric(x[1])]]
    comp_nodes_b <- Component_list[[as.numeric(x[2])]]
    sum(test_m[comp_nodes_a, comp_nodes_b])
  } else {
    NA
  }
})
preced_comp <- x_components_idx[which(x_comp_link > 0),]
colnames(preced_comp) <- c("precedent", "consequent")
ord_comp <- c(1:comp_count)                                                                                         # arrange components number in an array
ord_comp_UPDATE <- ord_comp
swap_track <- list()
for(i in 1:nrow(preced_comp)){
  if(i ==1){
    dummy_ord <- ord_comp
  } else {
    dummy_ord <- swap_track[[i - 1]]
  }
  a <- as.numeric(preced_comp[i,1])
  b <- as.numeric(preced_comp[i,2])
  a_pos <- as.numeric(match(a, dummy_ord))
  b_pos <- as.numeric(match(b, dummy_ord))
  if (a_pos > b_pos){                                                                                               # a should precede b but a is currently placed after b
    dummy_ord[b_pos] <- a
    dummy_ord[a_pos] <- b
  }
  swap_track[[i]] <- dummy_ord
  ord_comp_UPDATE <- dummy_ord
}
# obtain permutation
permuted_lines <- unlist(Component_list[ord_comp_UPDATE])                                                           # re-order components in the list
permuted_m <- test_m[permuted_lines, permuted_lines]
#### 03 - view the outputs ####
Component_list
permuted_lines
permuted_m

#### 04 - banchmark ####
library(igraph)
G_ig <- graph_from_adjacency_matrix(as.matrix(test_m), mode = "directed")
components(G_ig, mode = "strong")

test_SCC_tab <- lapply(1:length(Component_list),function(comp_count){
  x <- Component_list[[comp_count]]
  temp_rows <- length(x)
  Comp_ID <- rep(comp_count,temp_rows)
  cbind(Comp_ID, x)
})
tab_comp <- do.call(rbind, test_SCC_tab)
tab_comp_df <- as.data.frame(tab_comp)
#tab_comp_df[order(tab_comp_df$x),]
table(tab_comp_df$Comp_ID)