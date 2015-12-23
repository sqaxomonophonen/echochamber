# echochamber

A sound path tracer.

It renders an impulse response (suitable for convolution reverb) of a 3d model
by applying some of the principles used in a visual path tracer; rays are
traced from "microphones", bounced randomly off walls, in hope of hitting a
sound emitter. When a path to a sound emitter is found, a filtered impulse is
added to the impulse response where the filtering based on which materials were
encountered when bouncing off walls. The filtered impulse is delayed according
to distance traveled vs speed of sound.

I have no idea how "scientifically correct" any of this is, but it already produces convincing results, so uh.., yeah!

## Try this
```
./blender_export.sh blends/room2.blend room2.ecs
./ec init room2.ecs -I
./ec run room2.ecs
# ... wait a few seconds and press ctrl+c
./ec mixdown room2.ecs
mplayer room2.ecs.wav
```

# Blender integration

Cameras are microphones. Materials with non-zero emit parameters are sound emitters. Pretty much everything else
is currently ignored.


# Features
- Blender exporter
- Multiple microphones and emission groups (similar to light groups in LuxRender)
- "Indirect only" mode (`-I` switch to `ec init`); discards direct response, useful for dry/wet control

# TODO
 - Specular reflection, use a better BRDF function in general...
 - Sound emitter time offset (to allow multiple impulses at various times)
 - Various microphone characteristics (currently omni-directional only, how about cardioid and friends?)
 - Air scattering? (randomly filter/refract rays based on distance traveled)
 - Subsurface scattering?
 - Finish `ec mixdown`; one .wav per microphone; specify mix levels for each emission group
 - Ray-polygon intersection optimizations
 - Support for Blender smoothing groups?
 - Animation rendering? (render actual sound rather than just an impulse response)

# Thanks
 - Excellent compendium on global illumination: http://people.cs.kuleuven.be/~philip.dutre/GI/TotalCompendium.pdf

# License?
Nope, it's public domain.
