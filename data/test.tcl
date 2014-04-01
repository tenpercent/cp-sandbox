# create base cylinder 'c' with R=1, H=5
pcylinder c 1 5

# create sphere 's' with R=1.5
psphere s 1.5

# fuse cylinder and sphere: r = c + s
bfuse r c s

# create box 'b'
box b -2 -0.5 3 4 1 1

# cut box from previous result: r = r - b
bcut r r b

# display the result
donly r
fit

