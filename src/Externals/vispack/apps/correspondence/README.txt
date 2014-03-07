

linear_run - works.  Does 1-to-1 correspondences with full constraints.   Quite slow for big data.
linear_labels - works.  Does sparse correspondences with constraints that labels sum to one.   A little slow for big data.
linear_labels_resrict - tries to do a dense labeling and enforce unique labels through a penalty.  Does not work so well
linear_labels2 - works.  This is more efficient for big data - does not solve svd, but solves contraints directly.
linear_labels3 - Like linear_labels2, but more efficient matrix operations for bigger data
linear_labels4 - Like linear_labels3, but enforces the constraint that sum of labels to a point within an object < 1
linear_labels5 - Like linear_labels4, but threaded
