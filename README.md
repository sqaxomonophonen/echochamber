# echochamber

A sound path tracer.

It renders an impulse response (suitable for convolution reverb) of a 3d model by applying some of the principles
used in a visual path tracer; rays are traced from "microphones", bounced randomly off walls, in hope of hitting a
sound emitter. Each bounce pushes a filter onto a stack analogous to how colors are multiplied together in a visual
path tracer. When a path to a sound emitter is found, an impulse (i.e. `f[0]=1`, `f[x]=0` for `x!=0`) is filtered
through the filter stack, delayed according to distance traveled vs speed of sound, and finally added to the output
impulse response.

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
 - Material control from Blender
 - Replace biquads with FIRs?
 - Better filter control (currently too easy to pass biquad coefficients that amplify certain frequencies and cause the response to "blow up")
 - Sound emitter time offset (to allow multiple impulses at various times)
 - Various microphone characteristics (currently omni-directional only, how about cardioid and friends?)
 - Air scattering? (randomly filter/refract rays based on distance traveled)
 - Subsurface scattering?
 - Finish `ec mixdown`; one .wav per microphone; specify mix levels for each emission group
 - Ray-polygon intersection optimizations
 - Support for Blender smoothing groups?
 - Animation rendering? (render actual sound rather than just an impulse response)

# Thanks
 - Biquad coefficient calulator: http://www.earlevel.com/main/2013/10/13/biquad-calculator-v2/
 - Excellent compendium on global illumination: http://people.cs.kuleuven.be/~philip.dutre/GI/TotalCompendium.pdf

# License?
Nope, it's public domain.
