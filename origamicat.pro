;-----------------------------------------------------------------------------
; Copyright (c) 2016, Bridget L. Falck, Mark C. Neyrinck, & Nuala McCullagh
;
; Distributed under the terms of the Modified BSD License.
;
; The full license is in the file LICENSE, distributed with this software.
;-----------------------------------------------------------------------------

PRO origamicat, posfile,velfile,dir,fname,boxsize,mass,kpc=kpc

; Calculates halo properties for ORIGAMI halos and writes data to
; binary catalog file.

; INPUTS:
;
; - posfile: full path and file name of particle positions
; - velfile: full path and file name of particle velocities
; - dir: directory where .dat files stored
; - fname: label describing this run
; - boxsize: in same units as positions
; - mass: mass per particle in M_Sun/h

; KEYWORDS:
;
; - kpc: set this if position file is in units of kpc/h. Catalog is
;   written in units of Mpc/h.


IF keyword_set(kpc) THEN mpc = 1.e-3 ELSE mpc = 1

hidfile = dir + fname + 'hns.dat'
cenfile = dir + fname + 'cen.dat'
catfile = dir + fname + 'cat.dat'

npart = 0L
openr,1,posfile
readu,1,npart
xpos = fltarr(npart) & ypos = fltarr(npart) & zpos = fltarr(npart)
readu,1,xpos
readu,1,ypos
readu,1,zpos
close,1


numhalos = 0L
openr,1,hidfile
readu,1,numhalos
hmems = {np:lonarr(numhalos), pids:ptrarr(numhalos)}
nump = 0L
FOR i=0L,numhalos-1 DO BEGIN
    readu,1,nump
    hmems.np[i] = nump
    temp = lonarr(nump)
    readu,1,temp
    hmems.pids[i] = ptr_new(temp)
ENDFOR
close,1

print,'# halos = '+strn(numhalos)


openr,1,velfile
readu,1,npart
xvel = fltarr(npart) & yvel = fltarr(npart) & zvel = fltarr(npart)
readu,1,xvel
readu,1,yvel
readu,1,zvel
close,1

openr,1,cenfile
readu,1,numhalos
xc = fltarr(numhalos) & yc = fltarr(numhalos) & zc = fltarr(numhalos)
readu,1,xc
readu,1,yc
readu,1,zc
close,1

;stop

g = 4.45e-9 ; gravitational constant in km^2 Mpc / s^2 M_sun
rho200 = 1.e6/g ; M(<r)/Vol = 200 rho_crit -> M(<r)/r^3 = 100^3/g


;skip = 5   ; skip center particles
skip = 1    ; don't skip
volpart = boxsize^3./float(npart)
vcx = fltarr(numhalos)
vcy = vcx & vcz = vcx
m200 = fltarr(numhalos)
r200 = m200 & vmax = m200 & rmax = m200 & mmax = m200
maxr = fltarr(numhalos)
size = fltarr(numhalos)

print,'Calculating halo properties:'
FOR i=0L,numhalos-1 DO BEGIN

    ind = *hmems.pids[i]

    pos = transpose([[xpos[ind]-xc[i]],[ypos[ind]-yc[i]],[zpos[ind]-zc[i]]])
    pos += (pos LE -boxsize/2.)*boxsize
    pos -= (pos GT boxsize/2.)*boxsize

    r = reform(sqrt(pos[0,*]^2.+pos[1,*]^2.+pos[2,*]^2.))
    r = r*mpc ; get distance in Mpc/h!
    s = sort(r)
    s = s[skip:(hmems.np[i]-1)]

    v = (findgen(hmems.np[i]-skip)+skip)/r[s] ; v index is 0 to nmem-1-skip
;    rho = volpart*v*3./(4.*!pi*r[s]^2.)   ; ditto with rho
    rho = mass*v/r[s]^2. ; M(<r)/r^3
    max = where(v EQ max(v))
    max = max[0]
    vmax[i] = sqrt(g*mass*v[max]) ; sqrt(GM(<r)/r)
    rmax[i] = r[s[max]]
    mmax[i] = (max+skip)*mass

    two = where(rho GT rho200,count) ; index in s
    two = max(two)
    IF count EQ 0 THEN BEGIN
        r200[i] = 0.0
        m200[i] = 0.0
    ENDIF ELSE BEGIN
        r200[i] = r[s[two]]
        m200[i] = (two+skip)*mass
    ENDELSE

    vcx[i]=avg(xvel[ind]) & vcy[i]=avg(yvel[ind]) & vcz[i]=avg(zvel[ind])

    maxr[i] = max(r)

    sizexyz = max(pos,d=2)-min(pos,d=2)
    size[i] = max(abs(sizexyz))*mpc ; kpc to Mpc

ENDFOR

; Write catalog:

; get centers in Mpc
IF keyword_set(kpc) THEN BEGIN
    xc *= 1.e-3 & yc *= 1.e-3 & zc *= 1.e-3 
ENDIF

print, 'Writing catalog:'
openw,1,catfile
writeu,1,numhalos
writeu,1,xc
writeu,1,yc
writeu,1,zc
writeu,1,vcx
writeu,1,vcy
writeu,1,vcz
writeu,1,m200
writeu,1,r200
writeu,1,vmax
writeu,1,rmax
writeu,1,mmax
writeu,1,hmems.np
writeu,1,maxr
writeu,1,size
close,1 


END
