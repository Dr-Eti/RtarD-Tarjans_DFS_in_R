###  Find Strongly Connected Components in a digraph by Depth-First Search [Tarjan]

#### -1.0 - Revision notes ####
# Major edits on 30 10 2021 - 11 10 2021
#   A revision was required to investigate and fix a 'loss' of component elements. The loss emerged during the analysis of larger matrices than those provided in the examples below.
#   The main change introduced by this version is that the second lowlink update (the one during backtracking) 
#   may need few iterations before all nodes' lowlink values are actually set to correspond to the lowest lowlink value among their successors.
#   another, minor update is that, when the successor of a node is a sink, there is no reason for the node to inherit its successor's lowlink value, even if lower;
#   otherwise they might erroneously end up in the same component, which doesn't really make sense when the graph is directed.

#### -1.1 - intro: date, author, contact ####
# Coding timeline
#   Start:      13 07 2021
#   End:        06 08 2021
#   Fixes:      19 08 2021
#   for GitHub: 14 09 2021
#   This REV    02 11 2021
#
# Author:
#   Ettore Settanni
#   e.settanni@eng.cam.ac.uk
#   Cambridge
#

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
test_m <- matrix(c(
  0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0
),ncol = 13, byrow = TRUE)


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
c <- 0                                                                                                     # component counter

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
  
  # newly added
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
    sink_test <- which(is.na(w_labels))                                                                    # check if successor is a sink node
    # update predecessor (unless first node)
    # if (i > 1){                                                                                          # this old line was replaced: unless the first node is also a sink node, it may well be it is some other node's predecessor
    back_idx <- which(util_table_df[,"w_k_node"] == v)                                                     # ... I guess we could use this information to single out sink nodes already
    if (length(back_idx) > 0){  
      util_table_df[back_idx, "w_k_number"] <- node_numbering[v,"node_number"]
      util_table_df[back_idx, "w_k_lowlink"] <- node_numbering[v,"node_lowlink"]
    }
    if(length(sink_test) == 0){
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
          util_table_df[util_table_df_subset,][k,"w_k_node"] <- w_k                                        # could be removed... 
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
          util_table_df[util_table_df_subset,"w_k_node"] <- w_k                                            # could be removed... 
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
        test1 <- NA
        test2 <- NA
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
            dummy3 <- FALSE                                                                                # skip the next loop, break out and go increase k 
            dummy5 <- FALSE
          } else {
            dummy3 <- TRUE   
          }  
          # the test below is needed anyway unless we are in a sink node
          if(length(sink_test) == 0){                                                                      # we are not at a sink node
            # if successor is a sink then don't update or they will end up in the same component
            test3 <- FALSE
            successor_successors <- which(util_table_df[,"v_node"] == util_table_df[which(rownames(util_table_df) == v_label),"w_k_node"][k])
            if(util_table_df[successor_successors,"successors"][1]>0){test3 <- TRUE}
            if (test1 & test2 & test3){                                                                    # the current k successor has been visited PRIOR to v (has lower number), and it is on the stack already
              # update lowlink values - FIRST TYPE (going downwards)  
              node_numbering[v,"node_lowlink"] <- min(node_numbering[v,"node_lowlink"],node_numbering[w_k,"node_number"])
              util_table_df_subset <- which(rownames(util_table_df) == v_label)                            # filter for the current node v
              if(length(util_table_df_subset) > 1){
                util_table_df[util_table_df_subset,][1,"v_lowlink"] <- node_numbering[v,"node_lowlink"]    # always update the first line, corresponding to v
              } else {
                util_table_df[util_table_df_subset,"v_lowlink"] <- node_numbering[v,"node_lowlink"]
              }
              # update lowlink elsewhere
              #if (i > 1){
              other_idx <- which(util_table_df[,"w_k_node"] == v)
              if(length(other_idx)>0){
                util_table_df[ other_idx, "w_k_lowlink"] <- node_numbering[v,"node_lowlink"]
              }
            }
          }
          
          # Loop #4: revisit all the neighbours already visited
          while (dummy3) {
            if(n_successors !=0 & k !=0){                                                                  # we are not at a sink node / we have not re-visited all the neighbors
              # test if successor is a sink then don't update or they will end up in the same component
              test3 <- FALSE
              successor_successors <- which(util_table_df[,"v_node"] == util_table_df[which(rownames(util_table_df) == v_label),"w_k_node"][k])
              if(util_table_df[successor_successors,"successors"][1]>0){test3 <- TRUE}
              if (test2 & test3){                                                                          # the next node is already in the stack AND IT IS NOT A SINK
                # update lowlink values - SECOND TYPE (ascending)
                node_numbering[v,"node_lowlink"] <- min(node_numbering[v,"node_lowlink"],node_numbering[w_k,"node_lowlink"])
                util_table_df_subset <- which(rownames(util_table_df) == v_label)                          # filter for the current node v
                if(length(util_table_df_subset) > 1){
                  util_table_df[util_table_df_subset,][1,"v_lowlink"] <- node_numbering[v,"node_lowlink"]  # always update the first line, corresponding to v
                } else {
                  util_table_df[util_table_df_subset,"v_lowlink"] <- node_numbering[v,"node_lowlink"]
                }
                # update lowlink elsewhere
                #if (i > 1){
                other_idx <- which(util_table_df[,"w_k_node"] == v)
                if(length(other_idx)>0){
                  util_table_df[ other_idx, "w_k_lowlink"] <- node_numbering[v,"node_lowlink"]
                }
              }
              k <- k - 1                                                                                   # go back to the previous neighbour already visted
              if (k !=0){
                w_k <- w[k] 
                #re-do the test
                test0 <- is.na(node_numbering[w_k,"node_number"])                                          # we can jump onto this node (depth first)
                test1 <- node_numbering[w_k,"node_number"] < node_numbering[v,"node_number"]               # the next node has been visited already?
                test2 <- node_numbering[w_k,"node_onStack"] == 1                                           # the next node is in the stack already
              } else {
                test2 <- FALSE
                w_k <- (-1)
              }
            } else {                                                                                       # either we are at a sink node or we have re-visited all the neighbors
              dummy3 <- FALSE                                                                              # break loop 4 after moving one level up
              level <- level - 1                                                                           # go one back to the previous node in the depth-first sequence
              
              # The below wraps the second ("upwards") lowlink update into an iteration - without this, a few lowlink updates may go unnoticed
              if (level >= 0){
                
                # The below checks that, when you've reached level 0, all nodes' lowlink values matches the lowest  lowlink of their successors.
                # if that's not the case it may be that we need more rounds of update
                if (level ==0 & flag_extra_round){
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
                    level <- max(util_table_df[which(!is.na(util_table_df[,"level"])),"level"])           # start all over again, hoping we get all the lowlinks updated
                  } else {
                    flag_extra_round <- FALSE                                                             # break out of this
                  }
                }
                
                # GET BACK WHERE YOU LEFT THINGS AT THE PREVIOUS LEVEL
                # (the following fixes the case in which there is no successor going up one level)
                test_any_successor <- length(which(util_table_df[,"level"] == level & util_table_df[, "successors"] !=0))     # test if at this level there may be no successor
                while(test_any_successor == 0 & level > 1){
                  level <- level - 1                                                                                          # if there is no successor, keep going backwards
                  test_any_successor <- length(which(util_table_df[,"level"] == level & util_table_df[, "successors"] !=0))   
                }
                if (test_any_successor != 0 & level > 0){
                  back_idx <- which(util_table_df[,"level"] == level & util_table_df[, "successors"] !=0)                     # there may be more than one node with the same level if E.G. ONE NEIGHBOUR IS A SINK but the other isn't. We ignore the sink
                  v <- util_table_df[back_idx, "v_node"]                                                  
                  v_label <-  node_names[v]
                  n_successors <- util_table_df[back_idx, "successors"]
                  k_back_idx_a <- as.numeric(which(util_table_df[, "v_node"] == v))
                  util_table_df_subset <- util_table_df[k_back_idx_a,]
                  if(length(k_back_idx_a) > 1){
                    k_back_idx_b <- max(which(!is.na(util_table_df_subset[, "k_successor_idx"])))                             # retrieve last successor index explored for the current node
                    k <- util_table_df_subset[k_back_idx_b , "k_successor_idx"]
                    w_k <- util_table_df_subset[k_back_idx_b , "w_k_node"]
                  } else {
                    k <- as.numeric(util_table_df_subset["k_successor_idx"])
                    w_k <- as.numeric(util_table_df_subset["w_k_node"])
                  }
                  w <- as.numeric(util_table_df[k_back_idx_a, "w_k_node"])
                  w_labels <- node_names[w]
                  sink_test <- which(is.na(w_labels))                                                      # check if the successor is a sink node
                  # re-do test
                  test0 <- is.na(node_numbering[w_k,"node_number"])                                        # we can jump onto this node (depth first)
                  test1 <- node_numbering[w_k,"node_number"] < node_numbering[v,"node_number"]             # the next node has been visited already?
                  test2 <- node_numbering[w_k,"node_onStack"] == 1                                         # the next node is in the stack already
                }
              }
            }
          } # end of Loop 4 (dummy 3)
          if(level == 0) {
            dummy5 <- FALSE                                                                                # break loop 3
            dummy1 <- FALSE                                                                                # break loop 2 (ends up in the same place as if there were no successors)
          }
        }  # end of Loop 3 (dummy 5)
        #  end of "back tracking" sub
      } 
    } # end of Loop 2 (dummy1)
    
    # THIS IS WHERE WE POP THE STACK AND POPULATE THE COMPONENT
    # update count of nodes not numbered yet
    nodes_not_numbered_yet <- length(which(is.na(node_numbering[,"node_number"])))
    nodes_visited_so_far <- as.numeric(which(!is.na(node_numbering[,"node_number"]) & node_numbering[,"node_onStack"] == 1))
    if(!test0){
      Stack_S_BACKUP <- Stack_S
      Component_list_BACKUP <- Component_list
      # Start of "Pop nodes" sub, but only if we're done with depth-first search
      for(r in nodes_visited_so_far){
        x <- as.numeric(node_numbering[r,])
        if (x[1] == x[2]){                                                                                 # the node is a root node
          c <- c + 1                                                                                       # component number
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
          Component_list[[c]] <- to_delete
          # pop elements out of stack
          #idx_delete <- which(Stack_S == to_delete)
          temp_idx_delete <- rep(seq_along(Stack_S), sapply(Stack_S, length))                              # thread: https://stackoverflow.com/questions/11002391/fast-way-of-getting-index-of-match-in-list
          idx_delete <- temp_idx_delete[match(to_delete, unlist(Stack_S))]            
          Stack_S[idx_delete] <- NULL
          # update table
          node_numbering[to_delete,"node_onStack"] <- (-1)                                                 # mark elements that have been on stack, but no longer are
        }
      }
      # need to re-label the LEVELS used so that they do not interefere with the next component
      util_table_df_BACKUP <- util_table_df
      idx_change_level <- which(!is.na(util_table_df[,"level"]) & util_table_df[,"level"] > 0)
      util_table_df[idx_change_level,"level"] <- (-1)*util_table_df[idx_change_level,"level"]
      # restart for next component or terminate
      if(nodes_not_numbered_yet != 0){
        dummy0 <- FALSE                                                                                    # break Loop 1 (depths first exploration of successors, jumping between nodes) and move up one level
      }  else {
        dummy4 <- FALSE                                                                                    # break Loop 5: YOU'RE DONE
        dummy0 <- FALSE                                                                                    # break Loop 1
      }
    } ## end of "Pop nodes" sub
  } # end of Loop 1 (dummy0)
} # end of Loop 0 (dummy4)

#### 02.2 - Block triangular form (bonus) ####

# Based on Plain English description from: Hume D., Litsey J., and Plemmons (1981) p 272 - see this google book https://books.google.co.uk/books?id=pEMsAQAAIAAJ&pg=PA272&lpg=PA272&dq=duff+and+reid+block+triangularization&source=bl&ots=zlTT95Usx-&sig=ACfU3U3GkBU0av0mI459KObNAqyrc-1sTw&hl=en&sa=X&ved=2ahUKEwiWkKnBlKbyAhUSesAKHXmZArcQ6AF6BAgREAM#v=onepage&q=duff%20and%20reid%20block%20triangularization&f=false
# - all vertices within a strongly connected component are numbered sequentially
# - if there exists and edge from a vertex in component i to a vertex in component j, then all vertices in component i are labelled before all those in component j
n_components <- length(Component_list)
x_components_idx <- expand.grid(1:n_components, 1:n_components)       # equivalent to nested for
x_components_idx <- x_components_idx[, ncol(x_components_idx):1]      # flip
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
ord_comp <- c(1:c)                                                   # arrange components number in an array
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
  if (a_pos > b_pos){                                                # a should precede b but a is currently placed after b
    dummy_ord[b_pos] <- a
    dummy_ord[a_pos] <- b
  }
  swap_track[[i]] <- dummy_ord
  ord_comp_UPDATE <- dummy_ord
}
# obtain permutation
permuted_lines <- unlist(Component_list[ord_comp_UPDATE])                    # re-order components in the list
permuted_m <- test_m[permuted_lines, permuted_lines]
#### 03 - view the outputs ####
Component_list
permuted_lines
permuted_m

#### 04 - banchmark ####
library(igraph)
G_ig <- graph_from_adjacency_matrix(as.matrix(test_m), mode = "directed")
components(G_ig, mode = "strong")