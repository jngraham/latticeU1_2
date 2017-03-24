# latticeU1_2
An attempt to do the NPQFT problem set in C++.

We use a one-dimensional array with enough entries to hold an x link, a y link and a t link for each site. The variable "lattice" is a pointer to the zeroth entry in the array so we have no issues passing lattice to things because pointers are only 4 (or 8) bytes (not sure how many they take up with my compiler). If we passed lattice by refernce then we'd need to dereference twice to get a link, which is less than ideal. What we have is good.

We can think of the structure of the lattice as having blocks of entries corresponding to links in sites that all have the same time coordinate, and in each of those blocks we have blocks of entries corresponding to links in sites that all have the same y coordinate, with further blocks of entries for x coordinate. In each block of 3 entries corresponding to a site that has a certain (x,y,t), the zeroth entry is the link pointing away from the site in the x direction, the first entry the y link, and the second entry the t link. So, for example, we access the y link of a site at (x,y,t) by

```
lattice[3*x + 3*Lx*y + 3*Lx*Ly*t + 1]
```

 My file test.cpp is intended to make sure I know what I'm doing when it comes to indexing!

My file main.cpp updates the lattice the appropriate number of times and extracts our data from the configurations, writing them to .csv files at the end. From there I can analyze the data elsewhere (possibly with numpy and matplotlib... we'll see). The issue is making sure all the functions I call are correct, and that problem lies in operators.cpp and link_update.cpp.

Still not wholly confident that my update() function in link_update.cpp is correct, mainly because I rewrote it from scratch in a way to get rid of the huge amount of modular arithmetic I was doing and got a totally different answer for <cos(U_p)>, which is calculated using the avg_plaquette() function, which I did not change. I am confident that my avg_plaquette() is correct. But! I just rewrote avg_plaquette() in a way to stop evaluating % everywhere and I get a totally different <cos(U_p)>! I obviously have a long way to go...

I have yet to write the other operators. I'll get to that in due course.
