# MusicalHomology

## The Goal

I want to use tools from Topological Data Analysis to be able to cluster music by genre. The hypothesis is that songs have homological signatures that capture something akin to genre.

## Attempt #1

I use a discrete time fourier transform to decompose a .wav file into individual notes.

I then decompose the frequency into a polar decomposition where the radius is the octave and its angle is a function of the note. My logic here is that if I just compare frequency, I'll get that different C's are far apart from each other, when musically we would think of them as being similar.

Thus I now have a massive array of points in four dimensional space(2 for the note decomposition, 1 for amplitude, 1 for time).

I then attempted to run Persistent Homology on this data set(Which could be anywhere from 5k-50k points). This ended up being a computational nightmare.


## Attempt #2(Where I'm at now)

After a couple of months of studying Persistent Homology, I realized that it the algorithm amounts to a giant matrix reduction(nxn for n points). 

To try to reduce the amount of points I run K-Medioids to get each song down to 250 points. At this point the computations can finish but certainly not in any fast amount of time. The K-Medioids is the major computational bottleneck.

After using Persistent Homology on a number of songs to obtain Persistence Barcodes, I compute the Wasserstein metric between these barcodes.

The results are not what I was hoping for. A Paul Simon is marked significantly closer to a Mountain Goats song than it is to another Paul Simon song.


## What's next

A couple issues I see with my current implementation
- Persistent Homology is using Veitoris-Rips complexes to create the filtration, which is highly dependent on the distance between my points. Thus my implementation is dependent on the distances between Medioids, and not the actual clusters. I can try to compute some alternative metric to better capture the distance between clusters.
- The bottleneck distance may be a better fit for comparing Persistence Barcodes, but I need to do a deeper dive on the differences between the two metrics.
- If the above 2 solutions don't get me anywhere, I may try to implement the Mapper Algorithm to generate simplicial complexes from a song, and then just directly compute the homology of this complex.



## Have any ideas?

Please don't hesitate to [reach out](mailto:liorbd@outlook.com)
