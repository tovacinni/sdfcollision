# Mesh SDF Collsion

This is a small `libigl` based code snippet that performs mesh to SDF collision
using several different schemes, including the optimization-based scheme from
[Robust Optimization for Robust Signed Distance Function Collision](https://mmacklin.com/sdfcontact.pdf).

The implementation is naive in that it's CPU based. If I decide to one day use this code
beyond just as a learning exercise for myself, I'll probably implement this in CUDA. :)

