I gave it all I have. I hope this is helpful.

The Maxwell_2D_TE and Maxwell_Lorentz are more well doucmented.

If you run the err_script_exact.m, both methods should give you second order convergence.

Switching to pulse source is left for reader as an exercise.

The operator splitting folder is a little bit more chaotic. But overall structure is the same.

Avoid for loop like plague! Vectorize! Vetorize! Vetorize!
The speed of your code also dictates how many simulations you can run.
Faster code allows you to do your research faster! More experiments.
These methods should run on a laptop!
Remember! Grandpa and grandma done it with less RAM and CPU than you.

Develope your folder system. A lot of the code is repeated. So do it well, do it once.

With all my love and good wishes,
- Wei Xi 
