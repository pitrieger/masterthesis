#set a variable in R and save in a csv file

# Single One-liner
system('/Applications/MATLAB_R2021a.app/bin/matlab -nodisplay -r "a=2; b=1; display(a+b); exit"')

# Full MATLAB script
setwd("~/Desktop")
x <- 10
write.table(x, file='~/x.csv', sep=",", row.names=FALSE, col.names=FALSE)

#make a vector where each element is a line of MATLAB code
#matlab code reads in our variable x, creates two variables y and z, 
#and write z in a csv file
matlab.lines <- c(
  "x = csvread('~/x.csv')",
  "y=20",
  "z=x+y",
  "csvwrite('~/z.csv', z)")

#create a MATLAB script containing all the commands in matlab.lines
writeLines(matlab.lines, con="~/myscript.m")

#run our MATLAB script
system("/Applications/MATLAB_R2021a.app/bin/matlab -nodisplay -r \"run('~/myscript.m'); exit\"")
# read in the variable z we created in MATLAB
z <- read.table("~/z.csv")
z

#remove the temporary files we used to pass data between R and MATLAB
system("rm ~/x.csv ~/z.csv")

