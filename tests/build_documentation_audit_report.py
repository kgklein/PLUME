"""Assemble Stage 5A inventories and page-by-page read-only audit report.

Requires the artifacts from documentation_output_audit.py, isolated baseline
build/tests, documentation_workflow_probes.py and documentation_behavior_probes.py.
No pre-existing repository file is modified; only the requested audit report/JSON.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess
from documentation_output_audit import ROOT,BASE

# Semantics reviewed against functions.f90, vars.f90 and production call sites.
DEFINITIONS={
'heating_convention':('character(16)',"'corrected' | 'legacy'",'Diagnostic convention; corrected complete real-frequency energy and initialized same-frequency harmonics; legacy original calc_eigen arithmetic.','dimensionless label','Corrected collisionless -1/serial/2D heating; legacy all-kinetic and old collisional workflows.','Independent of Z/Bessel/normalization; legacy cold heating unavailable.'),
'z_backend':('character(16)',"'faddeeva' | 'legacy'",'Plasma dispersion Z backend, Z=i sqrt(pi) wofz; signed kpar convention retained.','dimensionless label','Every warm susceptibility; cold bypasses Z.','Independent of harmonics/diagnostic convention.'),
'normalization':('character(16)',"'rho_ref' | 'd_ref'",'Physical input length/temperature/velocity convention.','rho_ref uses perpendicular reference thermal gyroradius; d_ref uses inertial length','All modes; cold-containing requires d_ref.','Species 1 is physical reference; no separate reference selector.'),
'bessel_mode':('character(16)',"'adaptive' | 'fixed' | 'fixed128'",'Species-specific warm harmonic policy.','maximum absolute signed order N, sum -N..N','Warm collisionless/collisional dispersion; cold uses no Bessels.','fixed128 requires nfixed=128; FPC distribution uses separate fixed internal nbesmax=15.'),
'bessel_tol':('real','finite 0 < value < 1','Ordinary and n^2-weighted Bessel tail tolerance, not guaranteed tensor/root error.','dimensionless','Adaptive warm harmonics.','Both tails checked; arithmetic failure distinct from unresolved tail.'),
'bessel_nmax':('integer','1..46338','Adaptive maximum harmonic cap; warning and continue on unresolved tail.','absolute order','Adaptive warm harmonics.','Unresolved status recorded; cap warning once per cache state.'),
'bessel_nfixed':('integer','1..46338','Fixed maximum absolute harmonic order.','absolute order','Fixed warm harmonics.','fixed128 requires exactly 128.'),
'bessel_diagnostics':('logical','.true. | .false.','Enable .bessel cache/cutoff/convergence metadata.','logical','All calculations, no Bessels required for cold.','Legacy fixed128 cache counters are explicitly limited.'),
'betap':('real','positive for warm reference; exactly zero for cold reference','8*pi*nR*Tparallel,R/B0^2.','parallel reference beta','All modes.','d_ref requires tauS(1)=betap within 1e-8 relative; beta/reference-temperature scans synchronize them.'),
'kpar':('real','finite; nonzero in every wavevector scan and warm susceptibility','Parallel wavevector component.','k*rho_R or k*d_R','All calculations; signed drift/propagation.','All-warm kperp=0 and all-warm kpar=0 currently fail direct evaluation; see geometry issue.'),
'kperp':('real','finite; mixed/cold can be exactly zero','Perpendicular wavevector component, k=(kperp,0,kpar).','k*rho_R or k*d_R','All modes.','Cold/mixed exact parallel available; historical all-warm zero-kperp kernel singularity is not dispatched away.'),
'vtp':('real','positive','rho_ref: reference parallel thermal speed/c; d_ref: reference Alfven speed/c.','v_tparallel,R/c or vA,R/c','All modes; d_ref remains finite for cold reference.','No alias va_ref_c is an input; do not confuse internal auxiliary V with this physical a.'),
'nspec':('integer','must be positive; parser lacks explicit complete bounds check','Number of populations, in original indexed order.','count','All modes.','Provide every &species_j; species 1 fixes physical reference.'),
'nscan':('integer','nonnegative; 2 for option=2; option=4 one scan','Number of serial scans or 2D coordinate axes.','count','1/2/4; zero for supplied-guess single point.','Mixed/cold scans wavevector only.'),
'option':('integer','-1,0,1,2,4,5,6,7; 3 explicitly deprecated','Executable workflow selector.','enumeration','-1 direct; 0 root map; 1 supplied/map roots + serial scans; 2 grid; 4 repeated maps; 5 radial; 6 gyro FPC; 7 Cartesian FPC.','A zero process exit is not proof of valid outputs; deprecated/unsupported paths can STOP with status zero.'),
'nroot_max':('integer','positive; effective storage ceiling numroots=500, not fully input-validated','Number of supplied guesses / continued branches, not a guaranteed cap on .roots file lines.','count','Map searches, supplied roots, scans/FPC.','&guess_j required when use_map=false; map .roots writes all found refined candidates before some selection filters.'),
'use_map':('logical','.true. | .false.','Choose frequency map search versus supplied guesses.','logical','Options 1/2/5/6/7; option0 maps regardless.','Map needs &maps; guesses need indexed &guess_j.'),
'low_n':('logical','.true. | .false.','Select extra n=0,+1,-1 harmonic outputs in historical mode blocks.','logical','Heating scans and selected legacy susceptibility storage.','Corrected .energy always exposes all seven channels; standard sweep reader assumes low_n=true.'),
'collision_type':('integer','0 collisionless; 1 nonconservative Krook; other values STOP','Collision operator selection.','enumeration','All-kinetic only for finite rate.','Corrected heating requires 0. Legacy collisional diagnostics retained, with outstanding xz/yz questions; not new accuracy claims.'),
'kn':('real','positive physical mean free path ratio; incomplete finite/positivity validation','Neutral mean free path / selected reference length.','Kn_d=Kn_rho*sqrt(beta_R*alpha_R)','All-kinetic collisions.','Kn>1e25 gives nu=0; default1e30; nu_d=sqrt(mu*tau_d)/(sqrt(2)*Kn_d). nu_ns is derived, not an input.'),
'writeout':('logical','.true. | .false.','Controls much console verbosity, not a master output-file switch.','logical','All modes.','Python set_params currently incorrectly derives this from use_map rather than writeOut argument.'),
'dataname':('character(100)','existing relative data subdirectory','Subdirectory under data/ for output.','path','All output-producing workflows.','Fortran writers generally do not create it; mkdir -p data/<dataName> first.'),
'outputname':('character(100)','output basename','Common prefix for mode/maps/metadata outputs.','filename component','All modes.','Files may be overwritten; companion single/scan runs reset prior .energy.'),
 'taus':('real','rho_ref >0; d_ref >=0, exact zero cold','rho_ref: Tparallel,R/Tparallel,s; d_ref:2*Tparallel,s/(mR*vA^2).','inverse temperature ratio or Alfvenic temperature','Every species.','Exact zero is model selector only in d_ref; negative invalid; no small-temperature threshold.'),
'mus':('real','positive','Reference/species mass ratio mR/ms.','dimensionless','Every species.','Physical reference value should be1; parser warns only for mu1>1+1e-7, not symmetric/terminal.'),
'alphs':('real','warm >0; cold ignored','Species perpendicular/parallel temperature ratio.','dimensionless','Every species; cold placeholder irrelevant.','Physical alpha_s distinct from auxiliary alpha*=1; rho_R uses alpha_R.'),
'qs':('real','nonzero, signed','Reference/species charge ratio qR/qs.','dimensionless','Every species.','Gyrofrequency mu/Q signed; plasma prefactor uses Q^2. No neutrality hard constraint.'),
'ds':('real','nonnegative physical density; cold references/cold kernels enforce; warm-only validation incomplete','Species/reference density ratio ns/nR.','dimensionless','Every species.','Zero density has unavailable relative-density moments; not physically defined zero perturbation.'),
'vvs':('real','finite signed parallel drift','Equilibrium Us/vA,R along B0.','Alfven velocity','Every species, including species1 with warning.','No forced zero-current frame or automatic charge/current balancing.'),
'loggridw':('logical','.true. | .false.','Select logarithmic real-frequency map coordinates.','logical','&maps, option0/maps.','Confirmed endpoint-origin inconsistency; do not present intended log endpoints as current correct behavior.'),
'loggridg':('logical','.true. | .false.','Select logarithmic negative-imaginary map coordinates.','logical','&maps.','Code forces negative sign and uses incorrect origin with logarithmic increments; see P0-map.'),
'omi':('real','finite; >0 if log','Requested real frequency lower bound.','Omega_R','&maps.','Linear endpoints work; log origin currently exponentiates omi rather than log10(omi).'),
'omf':('real','finite; >0 if log','Requested real frequency upper bound.','Omega_R','&maps.','Bounds ordering/grid positivity not comprehensively validated.'),
'gami':('real','finite; nonzero magnitude if log','Requested imaginary frequency lower bound.','Omega_R, positive growth','&maps.','Log mode always negative; intended bounds currently not reproduced.'),
'gamf':('real','finite; nonzero magnitude if log','Requested imaginary frequency upper bound.','Omega_R','&maps.','Same log-origin issue as gami.'),
'nr':('integer','positive intended, not fully validated','Real-frequency map intervals, generating nr+1 nodes.','interval count','&maps.','Default128 ->129 nodes; grid has (nr+1)*(ni+1) rows.'),
'ni':('integer','positive intended, not fully validated','Imaginary-frequency map intervals, generating ni+1 nodes.','interval count','&maps.','Default128 ->129 nodes.'),
'positive_roots':('logical','.true. | .false.','Select candidates with real frequency > -1e-10 after map file writing.','logical','Map-selected continuation roots.','Does not mean positive growth; .roots written before final filter; hard gamma<1 filter also exists.'),
'g_om':('real','finite supplied frequency guess','Real part of supplied frequency.','Omega_R','&guess_j, use_map=false; direct option-1 uses first.','Secant roots are refined except direct evaluation.'),
'g_gam':('real','finite, signed','Imaginary part of supplied frequency; growth positive.','Omega_R','&guess_j.','Strongly damped/growing heating can be qualified/unavailable.'),
'scan_style':('integer','-1 compound k;0 globals;1..nspec individual species','Select family of scan coordinate.','enumeration','Options1/2/4; legacy parameter scans.','Mixed/cold allow wavevector scans only.'),
'scan_type':('integer','style-1:0 vector,1 angle,2 magnitude; style0:0 kp,1 kz,2 beta,3 vtp; species:0 tau,1 mu,2 alpha,3 Q,4 D,5 U','Select coordinate inside scan style.','enumeration','&scan_input_j.','Vector swi/swf are final kp/kz; angle uses swi as terminal angle; magnitude uses swf; not generic two endpoints.'),
'swi':('real','coordinate-dependent valid initial/target value','Scan range_i; compound vector/angle overload.','native coordinate; angle degrees','&scan_input_j.','For style-1 vector: target kperp; angle: target theta; initial physical point comes from current plasma state.'),
'swf':('real','coordinate-dependent valid final/target value','Scan range_f; compound vector/magnitude overload.','native coordinate; angle degrees','&scan_input_j.','For style-1 vector: target kpar; magnitude: terminal signed |k|.'),
'swlog':('logical','.true. | .false.','Logarithmic or linear continuation coordinate.','logical','Scans.','Log wavevector coordinates require positive endpoints; use linear for negative kpar.'),
'ns':('integer','positive for wavevector scans','Number of output intervals; normally ns+1 output points.','count','Scans.','Internal steps ns*nres; hidden substeps can fail without a written point.'),
'nres':('integer','positive for wavevector scans','Substeps per output interval.','count','Scans.','Writes when mod(internal_index,nres)=0; companion point is internal index.'),
'eigen':('logical','.true. | .false.','Request E/B and species density/velocity moments.','logical','Serial/2D; grid uses scan_input_2.','Maps currently accept but ignore this flag; Ex=0/degenerate gauge -> NaN.'),
'heating':('logical','.true. | .false.','Request species absorption and channels; corrected also companion energy.','logical','Serial/2D, grid scan_input_2.','Corrected collisionless wavevector only; map corrected request rejected, legacy map request ignored.'),
 'tensor':('logical','.true. | .false.','Request six independent complex susceptibility elements per species.','ordinary chi','Serial/2D output; grid requires BOTH tensor flags.','Map accepts and omits it; Python make_scan lacks a tensor argument but dict can supply it.'),
'nrad':('integer','positive required model-row count','Number of radial model rows.','count','Experimental option5.','Requires modelName0.rad and one species table per population.'),
'nrad_sub':('integer','positive intended','Radial interpolation substeps.','count','Option5.','Default1; interpolate in selected input variables, not generally invariant between normalizations.'),
'modelname':('character(100)','valid input model basename/path','Prefix of <modelName>0.rad globals and <modelName>s.rad species tables.','path','Option5.','Global rows radius,beta,vtp; species rows radius,tau,mu,alpha,Q,D,U.'),
'k_scan':('integer','0..6','Radial wavevector selection:0 fixed,1 kz,2 kp,3 magnitude,4 angle,5 kp/kz plane,6 angle/magnitude plane.','enumeration','Option5.','Selects different &k_range namelist layouts; cold-containing radial unsupported.'),
}


def parameters(data):
    booleans={'radial_heating','radial_eigen','rad_heat','rad_eigen','rad_log','rad_log_k','rad_log_par','rad_log_perp','rad_log_theta'}
    defaults={'heating':False,'eigen':False,'tensor':False,'nrad_sub':1,'elecdircontribution':0.}
    for p in data['input_parameters']:
        name=p['name'];groups={g['group'] for g in p['groups']}
        if name in DEFINITIONS:t,allowed,definition,units,modes,dependencies=DEFINITIONS[name]
        elif name in booleans:
            t='logical';allowed='.true. | .false.';units='logical';modes='Experimental radial option5.'
            definition=('Request historical radial heating.' if 'heat' in name else 'Request radial fields/moments.' if 'eigen' in name else 'Choose logarithmic interpolation for the named radial wavevector coordinate.')
            dependencies=('Corrected radial_heating rejected; use legacy all-kinetic; k_range rad_heat/rad_eigen are read but are not assigned to scan diagnostic fields (global radial flags control output).')
        elif name.startswith(('vperp','vpar','vx','vy','vz')):
            t='real';allowed='finite ordered velocity bounds';definition='Velocity-grid '+('minimum' if name.endswith('min') else 'maximum')+' for '+name[:-3]+'.';units='species/derived thermal velocity units, not native k-length convention';modes='FPC option6 gyro or7 Cartesian.';dependencies='All species kinetic; supply bounds explicitly; cold-containing FPC rejected. Header uses solver thermal values even for d_ref input.'
        elif name=='delv':t='real';allowed='positive';definition='FPC velocity-grid spacing.';units='thermal velocity';modes='FPC6/7.';dependencies='No comprehensive grid validation; zero spacing can cause invalid arithmetic.'
        elif name=='elecdircontribution':t='real';allowed='1 Ex only,2 Ey only,3 Ez only;0/other unchanged';definition='Debug selection of electric forcing components for FPC.';units='selector';modes='FPC6/7.';dependencies='Changes a debug diagnostic forcing, not dispersion. Not new mixed FPC support.'
        elif name in {'nk','nk2','kres','kres2','ntheta','thetares'}:
            t='integer';allowed='positive intended';definition='Radial k_range '+('resolution substeps' if 'res' in name else 'output interval count')+' for its designated coordinate.';units='count';modes='Option5, layout depends on k_scan.';dependencies='Defaults not assigned; nk/nk2 and kres/kres2 for Cartesian axes, ntheta/thetares for angle; model interpolation already separate.'
        elif name.startswith(('kperp_','kpar_','k_','theta_')):
            t='real';allowed='finite; positive if logarithmic';definition=('Initial/final angle' if name.startswith('theta') else 'Initial/final wavevector coordinate')+' in radial k_range layout.';units='degrees' if name.startswith('theta') else 'native selected k-length convention';modes='Option5; k_scan selects active fields.';dependencies='Suffix1 initial/fixed,2 terminal; other fields may be absent for that k_scan layout.'
        else:raise RuntimeError('Missing reviewed parameter definition '+name)
        p.update(data_type=t,allowed_values=allowed,physical_definition=definition,normalization=units,applicable_modes=modes,applicable_species='all kinetic only for collisions/FPC/radial; cold-containing collisionless d_ref otherwise',dependencies=dependencies,
            validation='Enumerations/positivity enforced where explicit in cited parser; no assumed blanket finite/bounds validation. See dependencies and prioritized input-validation finding.',
            spelling=next(x for g in data['namelist_groups'] for x in g['parameters'] if x.lower()==name))
        if groups & {'fpc'}:p['applicable_species']='All species kinetic, finite warm reference; cold-containing FPC rejected.'
        elif groups & {'radial_input','k_range'}:p['applicable_species']='All kinetic only; radial physics/continuation experimental.'
        elif groups & {'species'}:p['applicable_species']='Original population ordering; tauS decides kinetic/cold; d_ref required for any cold population.'
        else:p['applicable_species']='All supported profiles subject to the specific mode/dependency conditions above.'
        if name in {'normalization','z_backend','heating_convention'}:p['validation']='Exact case-sensitive enumerated string check at read_in_params; invalid value ERROR STOP.'
        elif name.startswith('bessel_'):p['validation']='validate_bessel_input: exact mode, finite0<tol<1, positiveinteger counts<=46338; fixed128 requires128; invalid ERROR STOP.'
        elif name in {'betap','vtp'}:p['validation']='Thermal adapter requires positive beta/vtp; cold reference requires exactlybeta0 andpositivea; d_ref tau1-beta equality1e-8relative. Explicit finite checks also occur in numerical kernels, not a blanket parser check.'
        elif name=='taus':p['validation']='Negative tau ERROR STOP; exactlyzero cold selected onlyd_ref; rho cold rejected. Reference consistency and collision/mode restrictions enforced in vars validators.'
        elif name=='alphs':p['validation']='Warm<=0 ERROR STOP; cold value deliberately not compared/read by thermal kernel. Finite validation can be deferred to guarded numerical routines.'
        elif name=='mus':p['validation']='<=0 ERROR STOP; reference mu>1+1e-7 printsnonterminalERROR; no symmetric selfratio enforcement.'
        elif name=='qs':p['validation']='Zero ERROR STOP; signed values retained. Neutrality/current onlywarn.'
        elif name=='ds':p['validation']='Negative cold/reference validation rejects; warm-only parser lacks nonnegative check; guarded response/root mayfail later.'
        elif name in {'kpar','kperp'}:p['validation']='Wavevector scan preflight requiresfinite andkpar!=0, including hiddensteps; mixed warmkpar0 explicitlyrejected; all-warm exactaxis reachesinvalidresponse ratherthan a universal earlygeometrycheck.'
        elif name in {'ns','nres'}:p['validation']='Wavevector scan positivecount andinteger-product overflow checks; older parameter/radial paths lesscomplete.'
        elif name=='nscan':p['validation']='2Drequires2,option4rejects>1 withSTOP; cold profilesrestrict modes; broadnonnegative/count initialization checks incomplete.'
        elif name in {'nspec','nroot_max','nr','ni','nrad','nk','nk2','ntheta','kres','kres2','thetares'}:p['validation']='No comprehensive earlycount/storage/order validation; validpositivevalues mustbe provided. numroots500 isinternalstorage, not validateduserupperbound.'
        elif name=='nrad_sub':p['validation']='Read default1; radial_scan silently clampsmax(1,nRad_sub).'
        elif name=='heating':p['validation']='Corrected rejects collisions,maps,plasma-parameter scans; legacyrejectscold-containingheatingscans; secondaxiscontrols2Dheat.'
        elif name=='tensor':p['validation']='No error for mapignored flag;2Drequiresbothaxisflags (source AND atopen/write),unlikeeigen/heatingsecondaxisrule.'
        elif name in {'loggridw','loggridg','omi','omf','gami','gamf','kn','delv'} or name.startswith(('vpar','vperp','vx','vy','vz')):p['validation']='Required validphysical values are not comprehensively checked at input; see log-map origin bug andFPC/Kn restrictions. Do not infer implementedvalidation fromintendedallowedrange.'
        elif name=='option':p['validation']='Executable SELECT CASE;3andunrecognizedvaluesSTOP andcanexit0; othermode restrictions invars/scanvalidators.'
        elif name in {'dataname','outputname'}:p['validation']='Fortrancharacter100 truncation; fileOPEN runtimeerrors ifdirectorymissing; noautomaticdirectorycreation.'
        else:p['validation']='Fortran namelist type parsing; no additional comprehensive range/finiteness/semantic check identified beyond the cited dependencies.'
        if name in defaults:p['default']=defaults[name]
        if not p['examples'] and groups&{'radial_input','k_range'}:p['working_test']='tests/test_normalization.py::test_radial_model_interpolation (synthetic .rad fixture)'
        if not p['examples'] and groups&{'fpc'}:p['working_test']='tests/test_normalization.py::test_fpc_diagnostics (gyro/cartesian normalized fixtures)'
    return data['input_parameters']


def capabilities():
    profiles=['kinetic_rho_ref','kinetic_d_ref','fully_cold_d_ref','mixed_warm_reference_d_ref','mixed_cold_reference_d_ref']
    rows=[]
    def add(mode,states,diagnostics,notes):
        for profile,status in zip(profiles,states):rows.append(dict(mode=mode,plasma=profile,status=status,diagnostics=diagnostics,conditions=notes))
    add('single supplied/direct oblique',['SUPPORTED']*5,'E/B,moments; corrected warm heating/channels and all-species energy when resolved. Fully cold kinetic heating NA.','-1 evaluates given w and requires a validated Ex=1 mode for fields;1,use_map=false,nscan0 refines guesses but writes no standalone mode file.')
    add('linear frequency map search',['SUPPORTED']*5,'determinant/root maps; no E/B/moment/heating files','option0 oruse_map=true; map .roots contains candidates beyond nroot_max; selected modes subject to additional filters.')
    add('serial wavevector scans',['SUPPORTED']*5,'E/B,moments, selected tensor output, corrected heat/channels/energy.','option1; nonzero kpar at every hidden and written point; negative kpar uses linear coordinates; branch crossings/poles can fail.')
    add('2D wavevector scans',['SUPPORTED']*5,'Same; second scan controls output; companion has outer/point context.','option2,2 axes; independent Cartesian kp/kz or ordered angle/magnitude; no perpendicular points.')
    add('repeated frequency maps',['PARTIALLY SUPPORTED']*5,'eigen/tensor accepted but silently omitted; corrected heating rejected; legacy map heating ignored.','option4,one scan. Use map/root schema only; log grid endpoint bug applies.')
    add('plasma-parameter scans',['PARTIALLY SUPPORTED','PARTIALLY SUPPORTED','UNSUPPORTED','UNSUPPORTED','UNSUPPORTED'],'Kinetic fields/moments and legacy heat; corrected parameter-scan heat rejected.','Style0 beta/vtp,style>0 species properties. d_ref beta/tau1 synchronized; simultaneous independent scans rejected.')
    add('radial scans',['PARTIALLY SUPPORTED','PARTIALLY SUPPORTED','UNSUPPORTED','UNSUPPORTED','UNSUPPORTED'],'Experimental kinetic legacy fields/heat; no corrected energy.','option5,model0.rad plus indexed species .rad; no public primary radial example; regression exercises a small interpolation case, not arbitrary solar-wind models.')
    add('exact parallel kperp=0',['UNSUPPORTED','UNSUPPORTED','SUPPORTED','SUPPORTED','SUPPORTED'],'Cold/mixed Ex=1 fields/moments and corrected heating/energy if applicable.','All-warm ordinary dispatch remains thermal singular at kperp=0 and direct probe returns NaN/status1. Cold/mixed uses the analytic parallel dispatcher; this distinction needs public documentation.')
    add('exact perpendicular kpar=0 single point',['UNSUPPORTED','UNSUPPORTED','PARTIALLY SUPPORTED','UNSUPPORTED','UNSUPPORTED'],'Fully cold X field when Ex=1; O mode can have Ex=0 and be unavailable.','Cold determinant O/X validated by prior independent tests. Warm kpar0 is unsupported; every wavevector scan forbids kpar0 including cold.')
    add('finite-rate collisions',['PARTIALLY SUPPORTED','PARTIALLY SUPPORTED','UNSUPPORTED','UNSUPPORTED','UNSUPPORTED'],'Legacy all-kinetic field/heating/FPC behavior retained; corrected heating not enabled.','collision_type1 nonconservative Krook. Finite-rate anisotropic/drifting xz/yz accuracy is NOT VALIDATED; zz independently corrected. No cold collisions.')
    add('FPC gyro/Cartesian',['PARTIALLY SUPPORTED','PARTIALLY SUPPORTED','UNSUPPORTED','UNSUPPORTED','UNSUPPORTED'],'Historical kinetic FPC/distributions/debug eigen/mom; no corrected mixed-FPC or adaptive FPC-distribution harmonic extension.','options6/7; finite-thermal reference; legacy normalization in headers; readers have hardcoded species labels/missing-value loss. Numeric equivalence tests do not validate every operator or velocity integral.')
    add('general finite-rate operator accuracy',['NOT VALIDATED','NOT VALIDATED','UNSUPPORTED','UNSUPPORTED','UNSUPPORTED'],'Do not infer validated heating/work from an executable collisional path.','Prior audit found anisotropic/drifting finite-rate xz/yz discrepancies; longitudinal zz independently corrected. No new collisional derivation here.')
    add('general FPC moment/energy closure',['NOT VALIDATED','NOT VALIDATED','UNSUPPORTED','UNSUPPORTED','UNSUPPORTED'],'Existing thermal FPC generation is available as above; arbitrary numerical grid/harmonic/closure accuracy is a separate scientific claim.','No mixed-FPC extension or new independent FPC-operator validation in this stage; fixed internal distribution harmonics differ from adaptive susceptibility harmonics.')
    return rows


def schemas():
    return [
 dict(id='stdout',pattern='stdout, option=-1',producer='disprels::test_disp',header='text labels',columns='not a single rectangular table',order=['w complex and displayed determinant complex','kperp,kpar','E(3),B(3)','density complex per species','velocity complex per species','Ps per species','first4 Ps_split only','corrected Htotal/HE/HB/Hspecies + statuses'],normalization='Native k; omega/OmegaR; determinant raw for cold reference, multiplied by solver_vtp^6 for warm reference; API determinant always raw.',failure='NaN fields with statuses; legacy direct defect conditions can exit0.',readers=['linfpc.test_disp uses file-based helper, not a general stdout schema'],notes='Corrected companion contains all7 channels; stdout not7.'),
 dict(id='map',pattern='data/<D>/dispersion_<O>[_<scan-tag>].map',producer='disprels::map_search via option0/use_map/4',header=False,columns=5,order=['omega_r','gamma','log10(abs(det Lambda))','Re det Lambda','Im det Lambda'],rows='(nr+1)*(ni+1),blank separator per real-frequency row',normalization='OmegaR; ordinary unscaled determinant (consistent within map)',failure='invalid samples NaN and masked from minima, not substituted physical finite tensors',readers=['np.loadtxt/gnuplot example_map.plt'],notes='Log endpoint bug; map root candidates not guarantee branch identification; maps may ignore eigen/tensor flags.'),
 dict(id='roots',pattern='data/<D>/dispersion_<O>[_<scan-tag>].roots',producer='map_search; NOT refine_guess alone',header=False,columns='6+6S',order=['kperp,kpar,beta,vtp,Re w,Im w','species1 tau,mu,alpha,Q,D,U; repeat original species order'],normalization='Native input convention, OmegaR',precision='es14.4',failure='failed refined candidates omitted',readers=['linfpc.compute_roots (returns only complex frequencies)','np.loadtxt'],notes='Actual map writes all found refined candidates before final selection; .roots15 rows in example_map_par despite nroot_max4. Published/local only-first-N claim incorrect.'),
 dict(id='mode_no_diagnostics',pattern='data/<D>/<O>_<parameter>_<range>.modeN',producer='om_scan/om_double_scan',header=False,columns='6+6S',order=['native kp,kz,beta,vtp,Re w,Im w','6S original species parameters'],precision='es15.6e3',failure='failed roots NaN frequency; parameters retained',readers=['load_plume_sweep(nspec=S,eigen=False,heating=False)','generic np.loadtxt'],notes='2D blocks separated by blank lines; generic reader flattens axes. Filename endpoints use scaled integer tags, can collide if values round to same tag.'),
 dict(id='mode_eigen',pattern='same .modeN with eigen=true',producer='scan diagnostic writer',header=False,columns='18+14S',order=['first6 native/root','Bxyz complex (6real)','Exyz complex (6real)','Ux,Uy,Uz complex per species (6Sreal)','dn complex per species (2Sreal)','6S parameters'],normalization='E=deltaE/deltaEx, Ex1; B=deltaB/deltaEx; U=deltaU/(c*deltaEx/B0); n=(delta n/n0)/(deltaEx/B0)',failure='unavailable Ex gauge -> fields/moments NaN; species failure isolates moments',readers=['load_plume_sweep explicit switches','SWIFT (reference velocity only, historical rho labels)','contours (schema inference can mislabel parameters as heat)']),
 dict(id='mode_heating',pattern='same .modeN,heating=true,low_n=true',producer='om_scan/om_double_scan',header=False,columns='6+14S without eigen;18+22S with eigen',order=['prefix6 or18+8S','Ps1..PsS','7 channels species1,then species2,...','6S parameters'],normalization='Ps=Im(E†chi_s(wr)E)/H;7 selected channels same H in corrected; legacy old H',harmonics=['n0yy','n0yz','n0zy','n0zz','fulln0','transverse+1','transverse-1'],failure='cold kinetic Ps/channels NaN; corrected failure/applicability NaN; fields/moments retained',readers=['load_plume_sweep(nspec=S,eigen=...,heating=True), current7stride fixed','load_plume_energy for statuses'],notes='0-based heat start18+8S ifeigen else6; total(s)=start+s; channel(s,c)=start+S+7s+c; parameters start=start+8S.'),
 dict(id='low_n_false',pattern='same .modeN,low_n=false',producer='legacy/corrected conditional writer',header=False,columns='eigen+heat18+15S+1;heat-only6+7S; legacy heat-only extra orphan Ew record',order=['usual prefix','Ps species totals','6S parameters','historical Ew final only eigen+heat, or extra legacy record'],normalization='Ew was unassigned, NOT the corrected H',failure='corrected final old Ew explicitlyNaN; legacy undefined/success can look finite',readers=['standard load_plume_sweep does NOT support this schema;4marker IndexErrors','explicit np.loadtxt with documented schema'],notes='Corrected orphan heat-only record fixed; legacy intentionally retained. No7 channels in old block, but corrected .energy retains all7.'),
 dict(id='susceptibility',pattern='data/<D>/<O>_tensor_<parameter>_<range>.modeN;2D<O>_tensor_<param1>_<param2>.modeN',producer='scan tensor writer/get_out_name',header=False,columns='4+18S',order=['kp,kz,beta,vtp','chi_xx for every species complex','chi_xy all','chi_xz all','chi_yy all','chi_yz all','chi_zz all','6S parameters species-major'],normalization='ordinary chi; no identity; transpose/sign symmetries supply other3',failure='NaN after invalid root',readers=['np.loadtxt explicit component-major schema'],notes='No frequency columns in this tensor file; associate with corresponding mode file/point. Not emitted by repeated maps even if requested.'),
 dict(id='energy',pattern='data/<D>/<O>.energy',producer='heating_diagnostics::write_energy_record, corrected -1/1/2',header='versioned plume-energy-v1; #columns and statuses/model metadata',columns='27+15S numeric + quoted source_file',order=['11int record/workflow/scan/outer/point/mode/root_valid/field_status/energy_status/applicability/response_code','16real kp,kz,wr,gamma,Htotal,HE,HB,Hspecies_total,fE,fB,fSpecies,fEabs,uncertainty,gamma_ratio,gamma_est,totalPs','15real per species Hs,fs,Qs(raw absorption),Ps,status,dmatrix,dproject,dstep,C1..C7','quoted source filename'],normalization='H=16piW/abs(deltaEx)^2; signed fractions unclamped; fEabs noncancelling contribution, not physical energy fraction',failure='explicit records for failed requested output roots; status8; cold heating status1 andNaN, energy mayvalid; strong/reactive flags preserve formalH butPsNA',readers=['load_plume_energy preserves quoted paths/failed rows'],notes='point internal integration index, not row when nres>1; outer/mode/source distinguish grid blocks; hidden substep failures are warnings only.'),
 dict(id='normalization_metadata',pattern='data/<D>/normalization_<O>.txt',producer='vars::write_normalization_metadata',header='key/value prose',columns='not rectangular',order=['normalization,units','composition/auxiliary scale/moment/diagnostic notes'],failure='stale d_ref removed on rho_ref name reuse',readers=['manual text; readers generally do not autojoin'],notes='Absent metadata does not reliably label arbitrary third-party/old data; input provenance still needed.'),
 dict(id='bessel_metadata',pattern='data/<D>/<O>.bessel',producer='adaptive_kinetic::finish_bessel_diagnostics',header='key/value lines plus#species header;no version marker',columns='10mixed fields per cachedspecies;optional6field unresolvedhistory',order=['species,lambda,N,capacity,converged,E0_bound,E2_bound,hits,misses,array_evaluations','optional unresolvedhistory:species,limit,lambda,tolerance,E0_bound,E2_bound'],failure='unresolved cap visible; only valid cached species have rows;cold/cache-unused species canbe absent;usecomposition metadata',readers=['np.loadtxt after header or explicit parser; no public primary wrapper'],notes='legacy fixed128 counters are not falsely reported as measured.'),
 dict(id='gyro_FPC',pattern='<O>.cpar/.cperp/.df1gyro.real/.df1gyro.imag.specieSS.modeMM',producer='fpc::compute_fpc_gyro',header='8labels+8values (tau,bi,kpar,kperp,vti,mu,wr*sqrtbeta/kpar,gi*sqrtbeta/kpar);5gridlabels+5gridvalues; separator',columns='matrix of velocity values; rows per vperp,columns per vpar',normalization='derived thermal solver parameters even for d_ref; correlate/distribution normalization per fpc source/README-JETPLUME',failure='old readers replaceNaN with0; no cold/mixed support',readers=['loadlinfpccepar,loadlinfpcceperp,loadlinfpcgyro_dist','linfpcplot gyro routines'],notes='specie02 labeled electron regardless physical species; no arbitrary-mixture interpretation; finite FPC harmonic nbesmax15 separate from adaptive susceptibility.'),
 dict(id='cartesian_FPC',pattern='<O>.cparcart/.cperp1/.cperp2/.dfs.real/.dfs.imag/.dfs3D.real/.dfs3D.imag.specieSS.modeMM',producer='fpc::compute_fpc_cart',header='labeled thermal parameters and grid bounds; projection blocks;3D files use separate header/axes',columns='multiple matrices/projection blocks,not a simple mode row',normalization='historical thermal/velocity units; sqrtbeta phase-speed header; fixedFPC physics',failure='old loadlinfpccart replacesNaNs with0 and recursively shifts header offsets on any exception',readers=['loadlinfpccart,loadlinfpccart_dist,loadlinfpc3d,loadlinfpc3d_dist','linfpcplot Cartesian routines'],notes='Header-offset recovery is not schema validation; 3D reader belongs to a separate explicit block schema.'),
 dict(id='FPC_debug',pattern='<O>.eigen.modeMM and<O>.mom.modeMM',producer='fpc routines, debugging',header=False,columns='eigen18+16S;mom11S',order=['eigen:6root/plasma,12fields,6Svel,2Sdens,S Ps,7Sselected','mom:2S density,6Svel,SjxEx,SjyEy,SjzEz'],normalization='historical calc_eigen/FPC debug conventions; not corrected .energy',failure='debug reader loadeigen hardcodes2species andPe=0; omitted drift density current/conjugation in derived products',readers=['loadeigen/loadmoms onlytwo-speciesdebug','do not use to establish general species energy closure'],notes='No species parameter suffix; thus NOT the same schema as .modeN scan despite similar prefix.'),
 dict(id='radial',pattern='data/<D>/<O>.modeN',producer='radial_scan/om_radial',header=False,columns='rootonly7+6S;eigen19+14S;heat8+7S;eigen+heat20+15S',order=['kp,kz,beta,vtp,complexw','optional fields/moments/Ps(totalonly)','6S parameters','old unassigned Ew when heating','radius LAST'],normalization='radius in solar radii; remaining native input/derived legacy conventions',failure='older continuation/error handling, no corrected metadata/status contract; NOT VALIDATED for arbitrary models',readers=['explicit custom parser; standard sweep reader is not radial-aware;nonheatingprefixcanparsebutdropsradius;heatingexpectswrong7channel schema'],notes='No public standalone primary radial example; tests cover synthetic interpolation only. Keep experimental.'),
    ]


def issues():
    data=[
('P0-01','P0','B','Scientific landing-page capability claim copied from ALPS','ford_project.md:44; public /','Serial nonrelativistic bi-Maxwellian/Krook and cold-fluid implementation; no MPI/OpenMP or relativistic kinetic equations in linked production modules. f192f70 introduced literal ALPS sentence/badge; later5bbca12 renamed subject.','Replace capability sentence with verified PLUME model; retain ALPS provenance in audit only.','small','None','low','yes'),
('P0-02','P0','A','Logarithmic map endpoint origins are inconsistent with increments','src/disprels.f90::map_search; behavior.json logarithmic-map','Requested wr .1..1 produces1.258925..12.58925; requested gamma -1e-4..-1e-2 produces -1.00023..-100.023.','Separately scoped map-coordinate correction with independent endpoint tests and explicit historical compatibility review; immediately warn in user docs.','small/medium','Define historical replay policy; no root algorithm change','medium: map roots change correctly','yes'),
('P0-03','P0','A','Exact-parallel all-kinetic dispatch is not supported despite broad propagation claims','disprels::disp/disp0 versus disp_mixed_parallel; behavior probes','All-warm kp=0 gives NaN determinant/status1, exit0. Cold/mixed use working analytic parallel path.','Document actual geometry matrix now; separately scope explicit warm parallel dispatch/validation; never add a kp floor.','medium','Independent parallel susceptibility/root benchmarks','medium','yes'),
('P0-04','P0','A','FPC readers convert unavailable numerical data to physical zeros','linfpclib/linfpc.py::loadlinfpccepar/loadlinfpccart','An injected NaN in an actual gyro FPC output is read as0; Cartesian source has same conversion.','Preserve NaNs/masks; make conversion an explicit legacy visualization option if retained. Do not claim mixed FPC support.','small/medium','Reader availability contract','low physics / medium plots','yes'),
('P0-05','P0','A','Resonance plotting reader invents heating columns in eigen-only files','contours/plot_plume_resonances.py::read_dispersion_output','Actual46-column warm-reference eigen-only file is marked heating_available; plotted first heating values .4 and1 are tau/mu inputs. Cold-reference file similarly reports0 and1.','Use input diagnostic flags and exact schema, reject ambiguous lengths; do not infer heating from residual species-parameter columns.','small','Schema inventory','low physics / medium plots','yes'),
('P0-06','P0','B','Published heating indices and scientific interpretation are obsolete','public page/output.html; current output.md differs','Published6stride and heating-only offset are wrong. Current active reader32marker configurations pass; public table still directs users to wrong species/channels.','Publish corrected7channel/offset schema; describe selected nonexhaustive channels and signed absorption, not irreversible heating.','small','Current schema authoritative','low','yes'),
('P0-07','P0','D','Markdown parsing corrupts normalization/scientific equations','public page/input.html normalization table; page/output.html species list; current local rendered build','Unescaped |Omega| splits table and discards d_ref expression; underscores become em elements (including ref_parallel attributes) inside raw dollar math. MathJax cannot restore stripped syntax.','Protect math with supported delimiters and escape table pipes; rendered HTML/visual equation tests, not just successful FORD exit.','medium','Pin renderer version','low code / medium docs','yes'),
('P1-01','P1','B','Published site and generated API lag current capabilities','112 published pages; only7modules vs current23generated modules','No published Faddeeva/Bessel/heating convention, cold/mixed fields/moments/energy guides;9current report pages404; source commit not exposed.','Publish versioned current user guide and API with visible source commit; development reports are archive, not user entry point.','medium','Correct build/link policy','low','yes'),
('P1-02','P1','D','Core navigation/download links broken','site.json 37broken link/anchor occurrences; current build64missing target occurrences','Tutorial installation/README, JET guide/notebook/logo, output-table SVGs, paired input downloads and report/result links404. Current workflow does not copy many linked assets.','Use generated-page-aware URLs; copy actual examples/assets; fail internal links/anchors in CI. External403 is access-limited, not proven broken.','medium','Page manifest','low','yes'),
('P1-03','P1','A','Wrapper loses valid Fortran controls on file roundtrip','plume_input::load_from_file; behavior probes','Drops low_n,loggridg,nr,ni,tensor,and all6Cartesian velocity bounds. Missing required loggridg/bounds can change or invalidate a later run.','Make schema-driven read/write roundtrip preserve recognized fields and reject unknown critical fields explicitly.','medium','Input inventory and roundtrip fixtures','medium wrapper compatibility','yes'),
('P1-04','P1','A','writeOut argument is ignored in Python constructor','plume_input::set_params;4combinations','writeOut is assigned from use_map:2of4combinations incorrect.','Assign requested argument independently; add narrow regression test.','small','None','low','yes'),
('P1-05','P1','B','Installation prerequisites and local doc instructions are inaccurate','INSTALL.md; Makefile; public install page','Build also needs C99 compiler and libm; ISO_C_BINDING/IEEE behavior beyond generic F90. Docs output is docs_out/index.html, not docs/index.html. https://https URL malformed.','Document tested GNU toolchain, bundled MIT Faddeeva, optional Python/plot/doc dependencies, existing data directory, actual entry points.','small','Toolchain inventory','low','yes'),
('P1-06','P1','D','Unpinned FORD builds succeed with ignored configuration and source parse warnings','workflow and FORD7.0.13 isolated log','srd_dir,base_url,ordered_subpage ignored; source comments rho_ref/d_ref parsed as unknown metadata; inactive gk_disp_gf parsing warning. Navigation order and some API units disappear.','Pin validated toolchain; fix metadata in a docs-focused task; exclude inactive source; validate expected API/page inventory and fail unexpected warnings.','medium','Renderer decision','low physics / medium docs','yes'),
('P1-07','P1','C','Historical introductory scans are poor clean-success smoke tests','example_map_par/high_beta_par; example_guess_double; cycle.sh','First two return0 with123/46root failures and many qualified/unresolved heating rows. Double grid timed out at20s with82derivative-failure records before termination. Tutorial default cycle now selects E plot, not documented map.','Retain as advanced stress examples with expected statuses; curate bounded supplied-root introductions and matching plot commands. Timeout is NOT proof of a broken solver.','medium','Independent branch benchmarks','low','yes'),
('P1-08','P1','B','Map .roots semantics differ from documented nroot_max cap','map_search writer beforefilters; actual example_map_par .roots','15refined candidate rows written despite nroot_max4. File is before final positive_roots/large-growth selection; supplied-guess single point writes no .roots.','Document file producer/selection/order precisely; distinguish candidate files, accepted/continued branches, stdout roots.','small','No code change required for accurate description','low','yes'),
('P1-09','P1','A','Map supplementary flags accepted but omitted','option4 with eigen/tensor/legacy heating','eigen/tensor requests produce only.map/.roots; corrected heating rejected; legacy map heating ignored.','Prefer clear validation/warning for ignored flags and document map-only support; preserve file formats.','small','Mode capability policy','low','yes'),
('P1-10','P1','B','low_n=false and historical Ew are not a general reader contract','disprels writers;4marker IndexErrors','Corrected writer fixes orphan heat-only Ew and sets old eigen+heat EwNaN. Legacy deliberately retains orphan/unassigned Ew. Active sweep reader only supports low_n=true.','State schema explicitly; add reader validation or a separate low_n=false parser; never advertise Ew as wave energy.','small/medium','Historical compatibility policy','medium for legacy parsing','yes'),
('P1-11','P1','B','Historical stage reports are mixed into current user documentation','normalization.md and many stage-prefixed input/output sections','Reports include earlier unavailability claims contradicted by later current implementation. Public normalization guide is explicitly a historical audit.','Keep archived reports with status/version banners; provide one current normalization/capability/diagnostic reference without chronology.','medium','Current definitions reviewed','low','yes'),
('P1-12','P1','A','Input validation incomplete for counts/reference ratios/nonfinite warm state','functions::spec_read/vars validation','No comprehensive counts/map-bounds/finite/warm density/Kn checks; species1 ratio mu>1 onlywarns, mu<1 notsymmetric, reference Q/D/tau_rho consistency notenforced.','Document required physical ratios/inputs; scope validation changes separately with legacy acceptance review, without mandatory neutrality/current constraints.','medium','Compatibility/design decision','medium','yes'),
('P1-13','P1','B','Three independent compatibility controls need one authoritative recipe','input.md earlyZ section,legacy_numerics.in,current heating section','Z+fixed128 reproduces dispersion numerics; heating also needs heating_convention=legacy. Old collisional tensor corrections are not reverted by numerical legacy flags.','Separate dispersion replay, heating replay, source-revision replay; publish all3controls plus revision when exact files required.','small','No physics change','low','yes'),
('P1-14','P1','B','Indexed namelist header grammar is stricter than Fortran namelist syntax','functions::input_unit_exist/get_indexed_namelist_unit; behavior probes','Exact lowercase &species_N/&scan_input_N/&guess_N must occupy their own stripped line; uppercase or trailing header comments are not found. Scratch files are written beside inputs.','Document accepted header grammar and writable input directory; separately improve parser validation/errors without changing physical parameter names.','small/medium','Input parser contract and compatibility tests','low/medium','yes'),
('P1-15','P1','B','JET guide coordinate prose and FPC amplitude/normalization contract need clarification','README-JETPLUME.md option6/7 prose; fpc::calc_fs1/compute_fpc_gyro/compute_fpc_cart','6 isgyro and7Cartesian in executable; guide prose interchanges descriptions. Gyro exbar=1,Cartesian exbar=sqrt(beta)/vtp; literal CE uses -Qs*v_i^2/2*Re(E_i*conj(df1/dv_i)). A dimensional power conversion is not independently established by the current normalization-equivalence tests.','Correct coordinate prose now; document literal thermal/raw amplitude conventions; defer full amplitude/charge/projection/velocity-integral validation to a dedicated FPC audit, without extending mixed FPC.','small docs; medium separate audit','FPC scientific definitions and independent velocity moments','low docs; unknown future numerical correction','yes'),
('P2-01','P2','A','Debug readers and duplicate legacy functions have misleading generality','loadeigen/loadmoms/load_plume_sweep_orig twice','loadeigen assumes2species andreturnsPe=0; derived j.E ignores conjugation/density drift; second original function shadows first. Debug sweep reader returns indices, not values.','Label debug-only; select one maintained generic reader; deprecate aliases with tests; do not treat fabricatedPe as physical heating.','medium','Reader policy','low physics / medium users','yes'),
('P2-02','P2','B','Plot utilities retain fixed species and rho labels','linfpcplot::plot_disp_power_2spec;gnuplot examples;SWIFT','Species1ion/species2electron hardcoded; axesrho even for d_ref; SWIFT k fields called *_rho andpositionsrho; valid generic numeric phase only if units interpreted correctly.','Document narrow contracts and numerical units; add normalization metadata selection later; no inference of new cold/FPC physics.','medium','Reader/plot schema contract','low','yes'),
('P2-03','P2','D','No numerical-example/link/schema CI despite installation claim','only doc.yml and ADS maintenance workflow','Docs runs onPR/push; no workflow invokes make/unittest; bibliography page omitted withoutsecret. Successful doc exit allows broken targets.','Add isolated build/test/example/schema/link checks and deterministic no-token bibliography placeholder. Do not run unbounded research scans in CI.','medium','Curated smoke set and pinned doc tooling','low','yes'),
('P2-04','P2','B','Run provenance and interrupted-run association need explanation','no headers on modes; .energy has root context not run UUID','4831completed companion records match mode rows, including failed roots. Timed-out run canleave flushed companion ahead of buffered mode rows. Unique output names prevent stale glob confusion.','Document companion joins/statuses and completion checks; consider later optional completion/run metadata without changing historic columns.','small/medium','Output schema','low','yes'),
    ]
    return [dict(id=x[0],priority=x[1],category=x[2],title=x[3],location=x[4],evidence=x[5],recommendation=x[6],scope=x[7],dependencies=x[8],regression_risk=x[9],independent=x[10]) for x in data]


def page_findings(site):
    result=[]
    for p in site['pages']:
        url=p['url'];path=url.removeprefix(BASE);title=p.get('title') or path
        if p['status']!=200:
            desc='Expected/current or linked page is not served.';behavior='HTTP '+str(p['status']);fix='Publish/copy canonical page or correct incoming link; do not call it reviewed content.';priority='P1';destination='Relevant user guide / Developer archive'
        elif path in ['', 'index.html']:
            desc='Landing page describes ALPS-like relativistic/parallelized model; duplicate Readme entry, unfinishedJETcitation.';behavior='Current linked PLUME executable is serial/nonrelativistic; cold/mixed models and newdiagnostics not introduced here.';fix='Use one concise model/quickstart landing page; remove inheritedALPS claims and placeholders.';priority='P0';destination='Home / Quick Start'
        elif path=='page/index.html':
            desc='Readme user overview predates mostcold/mixed features and new numerical/diagnostic defaults.';behavior='CurrentREADME source adds major capabilities; publicdownloads/JET links broken.';fix='Publish current model matrix and minimal first run; avoid second competing home.';priority='P1';destination='Overview / Installation / Quick Start'
        elif path=='page/input.html':
            desc='Legacy/dualnormalization input summary; table k normalization corrupt; lacksmodern controls.';behavior='89user-facing names/14namelistdefinitions;3independent diagnostic/numerical selectors; cold/mixed restrictions matter.';fix='Current input reference, explicitrequired/default values, geometriccapability notes, mathsafe tables.';priority='P0';destination='Input Reference / Normalization'
        elif path=='page/output.html':
            desc='Six-stride seven-channel indexing, outdated offsets/cap/scaling prose; brokenSVGtables.';behavior='Current7stride readers fixed; multiplemode/metadata/FPCschemas differ; .energy isseparate.';fix='Versioned exactschemas, availability/status definitions, readercontracts.';priority='P0';destination='Output Reference / Heating and Wave Energy'
        elif path=='page/tutorial.html':
            desc='Fourbranch wide map+scan and gnuplot walkthrough; installationlinkswrong.';behavior='Wide scans havefailed/qualified points, potentiallydifferentbranchselection; cycle default differs; notcleansuccess intro.';fix='One boundedverified supplied-root first run with explicitexpectedoutput, thenadvancedmapworkflow.';priority='P1';destination='Quick Start / Examples'
        elif path=='page/INSTALL.html':
            desc='Fortran-only dependencies and obsolete local documentationopenpath; malformedGitHubURL.';behavior='C99/ISO_C_BINDING, make,libm; runtimeoutputdirectories; currenttests201pass.';fix='Separate runtime/build/test/plot/doc requirements and testedtoolchains.';priority='P1';destination='Installation'
        elif path=='page/normalization.html':
            desc='Historical dualnormalization audit with olderzero-temperature restrictions.';behavior='Currentd_ref supportswarm/coldreferences; noauxphysicaltemperature; originalstageaudit isnotcurrentguide.';fix='Dedicatedcurrentnormalization definitions/table/pair; archivethisreport.';priority='P1';destination='Normalization / Developer archive'
        elif path=='page/citingpapers.html':
            desc='Scientificbibliography.';behavior='Notamodelcapabilityreference; buildwithoutADSsecret omitscorrespondingpage.';fix='Retain with datedprovenance anddeterministic no-secretfallback; bibliography completenessnotindependentlychecked.';priority='P2';destination='Citing PLUME'
        else:
            kind=path.split('/')[0];name=Path(path).stem
            desc=f'Generated {kind} reference for {name}; signatures/comments from older deployed source.'
            behavior='The name remains in source where applicable; deployed API lacks16current modules andnewexplicit-state interfaces. Historical comments do not guarantee current dispatch/units.'
            fix='Regenerate from pinned source/version; preserve useful signatures, classifyinternalutilities andverifyphysicalcomments.'
            priority='P2';destination='Developer / API Reference'
            if name=='zetout':
                behavior='Historical routine is now plasma_z_legacy::zetout_legacy; public zet_in selects default Faddeeva or explicit legacy.'
            if name in ['disp','disp0','disp1','calc_eigen','zet_in','zetout','read_in_params','spec_read','vars','specie','map_search','om_scan','om_double_scan']:
                priority='P1';fix='Regenerate and clarify backend/normalization/diagnostic routing; historical calc_eigen islegacy/FPC API, not corrected heating API.'
        sources=[]
        if path in ['', 'index.html']:sources=['ford_project.md']
        elif path=='page/index.html':sources=['README.md']
        elif path.startswith('page/'):
            candidate=Path(path).stem+'.md'
            sources=[candidate] if (ROOT/candidate).exists() else ['tests/'+candidate] if (ROOT/'tests'/candidate).exists() else ['workflow-generated or not staged']
        else:
            name=Path(path).stem.split('~')[0]
            for f in (ROOT/'src').glob('*.f90'):
                if re.search(r'(?im)^\s*(?:(?:pure\s+)?(?:complex\s+|real\s+|integer\s+|logical\s+)?)?(?:module|subroutine|function|interface|type)\s*(?:::)?\s*'+re.escape(name)+r'\b',f.read_text()):sources.append(str(f.relative_to(ROOT)))
            if not sources:sources=['FORD-generated index, historical name, or include-defined interface; use named module/procedure evidence']
        result.append(dict(url=url,source_path=sources,http_status=p['status'],title=title,existing_description=desc,verified_behavior=behavior,required_correction=fix,severity=priority,recommended_destination=destination,
            review='HTML text, structure, source mapping, links; no browser rendering' if p['status']==200 else 'Unavailable; HTTP status verified, content not reviewed'))
    return result


def curated(workflows):
    actual={Path(e['path']).name:e for e in workflows['examples']}
    rows=[
('01','Kinetic proton/electron and equivalent normalization pair',['normalization_rho_ref.in','normalization_d_ref.in'],'Retain physicalmatchedpair as one lesson;beta.4,alpha_ref1.5,κ(.05,.06),driftingelectrons.','Converted k_rho=k_d*sqrt(.4*1.5);vth/c=a*sqrt(.4);mode roots match retained precision.'),
('02','Low-frequency Alfven/kinetic-Alfven and sevenchannels',['heating_corrected_kinetic.in'],'Useexistingboundedrho_ref scan;identifybranchwithindependentpolarization.','Initial w≈.06098715236−1.2709968e-5i;H≈20342.6801;four n0channelsclose,7notexhaustive.'),
('03','Cold proton/electron parallel L/ion-cyclotron branch',['cold_proton_electron.in','cold_fields_parallel.in'],'Retainparallelroot/eigen companion;interpretR/Lwithsignedcharge.','Finite-electron-mass cold circular dispersion;w≈.06180521 atκpar.06 forfieldexample.'),
('04','Fullycold exactperpendicular X branch',['cold_proton_electron_perpendicular.in'],'Keepassinglepoint;neverclaimperpendicularscans.','IndependentO/Xdeterminant;warnthatOmodeEx0fieldsunavailable.'),
('05','Warmions/coldelectrons with correctedheating',['heating_corrected_warm_reference.in'],'Retain boundedwavevector scan.','Initialw≈.06153841958−5.5643255e-6i;H≈20384.8001;electronPs/channelsNaN,energyfinite.'),
('06','Coldreferences/warmdriftingelectrons',['heating_corrected_cold_reference.in'],'Retainboundedcase;adddocumenteddriftvariationonlyinfutureexampletask.','Initialw≈.05916734912−1.5887687e-6i;H≈20516.8792;betap0;auxscaleisnotphysicaltemperature.'),
('07','Interleavedmultispecies parallel drifts',['mixed_moments_multispecies.in'],'Retain moment/current/continuity lesson;orderingoriginal,kineticvsfluid.','FiniteE/B,momentsforall4species;species-currentmustincludedensity-driftterm;usesDs1.2/currentbalance.'),
('08','Anisotropydriven weak instability',[],'Proposednewinputfrom existingvalidated test_heating_energy_audit BENCHMARKS weak_growth;notcreatedinthisaudit.','τp.8,alpha_p2,τe.2,alpha_e1,κpar.3,κperp.05,a.01;w≈.3228628391+.00167495808i.'),
('09','Finite-mass whistler/Rmode',[],'Proposednewbounded coldorwarmrootbenchmarkfrom fullycoldtests;mustvalidate branch before publishinginput.','RecoverRmode n²=R;electronsfiniteμ1836;showpolarization/handedness andmatchingindependent root.'),
('10','Small 2Dkinetic/mixed wavevector grid',[],'Proposedboundedversion ofexisting 2Dtestfixtures,notlarge example_guess_double;skipknownpoles/crossings.','4–9writtenpoints withsingle-point roots/Efields/heating comparison;nonzerokpar;failedrecordsdemonstratedexplicitly.'),
('11','Vacuum electricenergyimportance',['heating_corrected_cold.in'],'Keepfullycoldenergyinput asordinarycase;proposecontrolledmoderate-a massratio4benchmark variant fromtest_electric_vacuum_substantial.','Moderate-a benchmark a.3,κpar.2,κperp.1,wr7.461100673,H53675.3078,fE≈.500000281;artificialmassratio mustbeexplicit.'),
('12','Exacthistoricalreproduction',['heating_legacy.in'],'Retainall3independentflags+sourcecommit;keepadaptive/Faddeevastressinputsassecondarynumericslessons.','heating_conventionlegacy,z_backendlegacy,bessel_modefixed128;exactpinnednumericalregressions.'),
    ]
    return [dict(slot=x[0],topic=x[1],existing_files=x[2],recommendation=x[3],benchmark=x[4],execution=[dict(file=f,returncode=actual[f]['returncode'],timeout=actual[f]['timeout']) for f in x[2] if f in actual]) for x in rows]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--artifacts',type=Path,required=True);args=parser.parse_args();p=args.artifacts
    data=json.loads((p/'static.json').read_text());site=json.loads((p/'site.json').read_text());workflow=json.loads((p/'workflows.json').read_text());behavior=json.loads((p/'behavior.json').read_text());docbuild=json.loads((p/'docbuild.json').read_text())
    inventory=parameters(data);findings=issues();pages=page_findings(site);matrix=capabilities();outputs=schemas();collection=curated(workflow)
    initial=json.loads((p/'before.json').read_text());changed=[name for name,h in initial.items() if not (ROOT/name).exists() or hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=h]
    baseline=(p/'tests.log').read_text();audit=(p/'audit-tests.log').read_text()
    count=re.search(r'Ran (\d+) tests in ([\d.]+)s',baseline);newcount=re.search(r'Ran (\d+) tests in ([\d.]+)s',audit)
    runs={r['path']:r for r in workflow['examples']}
    for e in data['examples']:
        e['execution']=runs.get(e['path'],dict(status='not run'))
        e['new_user_suitability']='secondary/advanced' if e['path'].split('/')[-1] in ['example_map_par.in','example_map_high_beta_par.in','example_map_perp.in','example_guess_double.in','adaptive_bessel.in'] else 'candidate bounded example; interpret statuses and branch identity'
        e['intended_problem']=(ROOT/e['path']).read_text().split('&params')[0].strip() or 'Historical example; physical branch is not proved by its filename.'
    source_docs=[]
    for f in sorted(ROOT.glob('*.md')):
        text=f.read_text();source_docs.append(dict(path=f.name,bytes=len(text.encode()),headings=re.findall(r'(?m)^#{1,3} (.+)$',text),published_canonical='page/index.html' if f.name=='README.md' else '/' if f.name=='ford_project.md' else 'page/'+f.stem+'.html',
          purpose='historical implementation report, not current user reference' if f.name=='normalization.md' else 'current/partially chronological user documentation',
          issues=['See prioritized findings and page-specific review.']))
    roadmap=[dict(stage='5B.0 targeted correctness contracts',scope='Independently fix/log-test map coordinate origins; explicit geometry validation/parallel support decision; retainNaNs inreaders; fix contour falseheating; preserve legacy physics/goldens.',acceptance='Independent grid endpoints; validparallel reference orclear rejection; unique markers; missingvalues preserved; unchanged prior numericalgoldens.',depends='No documentation restructuring required.'),
      dict(stage='5B.1 authoritative current user reference',scope='Model/capability overview;installation/quickstart;89input definitions;normalizationpair;diagnostics andlimitations;three-control legacy recipe.',acceptance='Every supported option/flag has a source/test reference; all primarycommands generate documented outputs with statuses.',depends='Document actualbehavior immediately; claim fixes only after5B.0.'),
      dict(stage='5B.2 public output/reader workflow',scope='Exactschema tables/typedreaders, explicitflags,low_n=false handling,wrapperroundtripandwriteOut;debug-onlyutilitylabels.',acceptance='32markerconfigurations,all speciesparameters/moments/channels,quotedenergy paths,NaN/failedroot/multimode/grid joins; no numericcolumn reordering.',depends='Schema inventory andcompatibilityreview; independent frompage styling.'),
      dict(stage='5B.3 documentation build/navigation',scope='PinFORD/renderer;properconfiguration/page metadata;protectmath/pipes;copyassets/examples;archiveoldreports;includeheatingguide;visiblecommit;linkchecker.',acceptance='BuildwithoutADSsecret;no unexpectedwarnings;alluserlinks/anchors/assets valid;renderedmath table tests andmanualbrowserinspection.',depends='Page manifest;canproceed alongside readercorrections.'),
      dict(stage='5B.4 curated examples andCI',scope='12lessoncollection fromthisaudit;boundedbranchbenchmarks;optionalnumericsstress/regression;safeoutputdirectories;doc/example/schema/makeunittestchecks.',acceptance='Freshcheckout first-run works;targetedexpectedroot/phase/statuschecks;timelimits distinguishunfinishedfromfailed;allcoretests/legacyregressions pass.',depends='Verifiedbranchbenchmarks and5B.1/5B.3;nomixedFPCclaims.')]
    created=['tests/documentation_output_audit.py','tests/documentation_workflow_probes.py','tests/documentation_behavior_probes.py','tests/build_documentation_audit_report.py','tests/test_documentation_output_audit.py','tests/PLUME_DOCUMENTATION_OUTPUT_AUDIT.md','tests/PLUME_DOCUMENTATION_OUTPUT_AUDIT.json','tests/PLUME_DOCUMENTATION_SITE_INVENTORY.json']
    summary=dict(stage='5A read-only documentation/input/output/workflow audit',utc=datetime.now(timezone.utc).isoformat(),commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
      initial_tracked_changes=False,initial_untracked_summary=dict(files=sum(x.startswith('??') for x in (p/'worktree-status.txt').read_text().splitlines()),note='Unrelated research artifact paths/contents omitted from public inventory; full start snapshot retained only in /tmp.'),protected_files=len(initial),changed_preexisting_files=changed,
      baseline=dict(tests=int(count[1]),seconds=float(count[2]),passed=bool(re.search(r'^OK\s*$',baseline,re.M)),optional_mpmath=True,location='isolated /tmp checkout with read-only .git link'),
      audit_tests=dict(tests=int(newcount[1]),seconds=float(newcount[2]),passed=bool(re.search(r'^OK\s*$',audit,re.M))),
      website=dict(base=BASE,retrieved_utc=site['retrieved_utc'],successful_html=sum(q['status']==200 for q in site['pages']),probed_html=len(site['pages']),missing_html=sum(q['status']==404 for q in site['pages']),broken_internal_link_occurrences=len(site['broken_internal_links']),
       rendered_browser=False,source_commit_exposed=False,last_modified_home=next(q['last_modified'] for q in site['pages'] if q['url']==BASE),inventory='PLUME_DOCUMENTATION_SITE_INVENTORY.json'),
      build=dict(make='success in temporary checkout; ordinary C99 compiler + gfortran + make + libm; no runtime third-party install',ford='7.0.13, requests2.31.0, graphviz; isolatedvenv;exit0withwarnings',generated_html=len(docbuild['pages']),missing_target_occurrences=len(docbuild['broken_internal_paths']),warnings=docbuild['warnings'],
        workflow='.github/workflows/doc.yml, push/PRmain;deploypushonlygh-pages;no numericaltestworkflow',limitations='No browser/JS renderer or ADS token used; published deploymentcause/commit notinferred.'),
      provenance=dict(landing_text='f192f70 initially copied an explicit ALPS sentence/buildbadge;5bbca12 laterchangedsubject toPLUME withoutremovingrelativistic/parallelizedclaims.',git_introduction='f192f70d245d6bdd97b4240a3f03efdd4dcff683',alps_reference='https://danielver02.github.io/ALPS/'),
      input_parameters=inventory,namelist_definitions=data['namelist_groups'],capability_matrix=matrix,documentation_pages=pages,repository_documentation=source_docs,output_schemas=outputs,
      scientific_uncertainties=dict(FPC=dict(literal_correlation='-Qs*v_i^2/2*Re(E_i*conj(df1/dv_i)); thermal species velocities',f0hat='exp[-(vpar-hatU_s)^2-vperp^2/alpha_s]',hatU_s='u_s*sqrt(tau_rho_s/(mu_s*beta_ref))',gyro_exbar='1',cartesian_exbar='sqrt(beta_ref)/vtp_solver',fs1_prefactor='-i*sqrt(mu_s*tau_rho_s/beta_ref)*(exbar/Qs)*f0hat/vtp_solver',projection='Cartesian sums omitted coordinate timesdelv;gyro integratesphi withvperp Jacobian anddelv/4centralcell substitution',physical_power_conversion='NOT VALIDATED in this documentation stage; no common gyro/Cartesian amplitude or dimensional calibration asserted',phase_header='literal w*sqrt(beta_ref)/Kpar; physical Alfven phase speed adds sqrt(alpha_ref) for historical perpendicular rho')) ,
      reader_compatibility=workflow['readers'],python_and_plotting=data['python_and_plotting'],behavior_probes=behavior,examples=data['examples'],
      regression_input_fixtures=[dict(path=str(f.relative_to(ROOT)),purpose='Pinned/regressionfixture,notprimaryuserexample',execution='coveredbybaselinewherecalled') for f in sorted((ROOT/'tests').glob('*.in'))],
      curated_examples=collection,prioritized_issues=findings,roadmap=roadmap,created_files=created,
      limitations=['All146probedHTMLpagesareinventoried;112successfulpagesexaminedasHTML/text,withsource-mappingforAPIanddetaileduser-guide/scientificentry-pointreview;noJavaScriptorvisualbrowserrendering.',
       'Externalthirdpartylinksnotexhaustivelyfetched;403/unreachablecannotbeclassifiedasdeadfromHTMLalone.','25of26projectexamplescompletedwithin20s;oneheavygridinterrupted:timeoutnotproofinvalidphysics.',
       'Noindependentbranchidentificationforallhistoricalwide-scanroots;manyqualifications/failuresrecorded.','Untrackedresearchinput/plot/outputtreesnotpublishedorinterpreted;all26inputs/examplefilesandtrackedutilitiesincluded.',
       'Finite-rateoff-diagonalcollisionaccuracyandnewFPCphysicswerenotre-derived;priorauditlimitationsretained.'])
    (ROOT/'tests/PLUME_DOCUMENTATION_OUTPUT_AUDIT.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    (ROOT/'tests/PLUME_DOCUMENTATION_SITE_INVENTORY.json').write_text(json.dumps(dict(published=site,current_local_build=docbuild),indent=2)+'\n')
    render(summary,ROOT/'tests/PLUME_DOCUMENTATION_OUTPUT_AUDIT.md')
    if changed:raise RuntimeError('Pre-existing files changed: '+str(changed))
    print('Audit assembled; protected pre-existing files unchanged:',len(initial))


def table(headers,rows):
    def cell(x):return str(x).replace('|','\\|').replace('\n',' ')
    return '| '+' | '.join(headers)+' |\n|'+'|'.join(['---']*len(headers))+'|\n'+''.join('| '+' | '.join(cell(x) for x in row)+' |\n' for row in rows)+'\n'


def render(d,path):
    issue=d['prioritized_issues'];site=d['website'];old=d['baseline'];new=d['audit_tests']
    text=f'''# PLUME Stage 5A: documentation, input, output and workflow audit

## 1. Executive summary

**Read-only audit. No pre-existing production source, reader, example, documentation,
plotting script, build/deployment configuration or test expectation was changed.**
The numerical work is committed at `{d['commit']}`. Tracked files were clean at
audit start; existing untracked research artifacts were preserved. The audit
checks {d['protected_files']} pre-existing tracked/source/test/example hashes;
all remain unchanged. Only the eight audit artifacts listed below were added.

The solver's current public contracts are substantially ahead of the deployed
user guide. Corrected collisionless kinetic/mixed heating, complete wave energy,
cold-reference mixtures, physical species moments, Faddeeva Z and adaptive
harmonics exist, but the published user pages mostly describe the earlier
normalization stage. A successful build or test suite does not establish the
correctness of every workflow or every scientific description.

The most significant **confirmed** findings are:

- The PLUME landing-page relativistic/parallelized sentence was carried over
  from ALPS. This is not confusion between websites: commit `f192f70` introduced
  the sentence with the subject ALPS and its build badge; `5bbca12` later renamed
  the subject. Current linked PLUME equations are nonrelativistic and the executable
  has no MPI/OpenMP path. Correct the PLUME description without attributing ALPS
  capabilities to it.
- Logarithmic map coordinates do not reproduce their requested bounds. A 3x3
  probe requesting wr=.1..1 produces **1.258925..12.58925**; gamma=-1e-4..-1e-2
  produces **-1.00023..-100.023**. The increment uses log10 bounds but the origin
  does not. This is a separate code correction, not an input-doc wording fix.
- Exact-parallel **all-kinetic** dispatch remains singular at kperp=0, returning
  NaN/status1 with process exit0. Cold/mixed exact-parallel dispatch works.
  Do not advertise a universal parallel limit based on tests of the mixed path.
- The resonance plotter reads species parameters as heating in eigen-only files:
  actual false “heating” values include **tau=.4 and mu=1**. Old FPC readers also
  turn NaN into **zero**. These are interpretation defects, not a plotting-style issue.
- The published seven-channel output indices are stale, while the **current active
  sweep reader is already fixed**. All 32 unique-marker current/debug combinations
  pass, with the debug reader correctly interpreted as returning indices. Cold/failed
  NaNs and **4,831** completed energy-companion records retain correct associations.

The isolated existing suite passes **{old['tests']} tests** in {old['seconds']} s,
including exact legacy regressions and checked builds. Eight additional frozen
characterization tests pass ({new['seconds']} s). No golden values were edited.
All 26 public input examples were attempted in fresh temporary run directories:
25 completed within the 20 s limit; the large 2D example was interrupted. Several
wide historical scans complete with failed roots and inapplicable/unresolved
heating, so “exit0” is not a scientifically sufficient success criterion.

## 2. Repository, build and published-site status

The executable is `plume.e`, built by the ordinary Makefile using GNU make,
Fortran with the repository's promoted double-precision/IEEE configuration,
an ordinary **C99 compiler**, and libm. Johnson's Faddeeva code is bundled;
there is no Python/SciPy/runtime-library requirement for the executable.
Test/reader use adds Python, NumPy, SciPy and Matplotlib; optional mpmath enables
arbitrary-precision references. GNU Fortran11.4 was exercised here; other
compiler flags are not automatically validated merely by being listed in Makefile.
The repository license is BSD-2-Clause; bundled Faddeeva has its own MIT notices.

`./plume.e path/to/file.in` is the entry point. The parser expects the `.in`
suffix and has a 100-character argument buffer; no robust help/invalid-suffix
workflow exists. Indexed species/guess/scan headers must be exact lowercase
names on their own line: uppercase or trailing header comments fail the manual
extractor despite ordinary Fortran field names being case-insensitive. The parser
writes a `.scratch` file beside the input, so that directory must be writable.
Create `data/<dataName>` before running: Fortran generally
opens output paths but does not create them. `writeOut` is verbosity, not a
master switch preventing files. Supplied-root option1/nscan0 prints roots,
but it does **not** itself write a standalone `.roots` file. Map search does.

The published site was accessible. We fetched/inventoried **{site['probed_html']}
HTML URLs**, of which **{site['successful_html']} returned200** and **{site['missing_html']}
returned404**. Counts include deliberately probed current-but-unpublished guide/API
URLs, not just navigation targets. The reachable published graph has112successful
HTML pages (including aliases, lists and API pages); each has a page-by-page record
below. There are **{site['broken_internal_link_occurrences']} broken internal link/anchor
occurrences**, including repeated occurrences of the same target. The homepage
Last-Modified header was `{site['last_modified_home']}`. The site does not expose
an unambiguous source commit; no deployment-failure cause or exact deployed
revision is inferred from staleness alone.

Every reachable user-guide page was reviewed as text/HTML. Generated API pages
were examined through complete HTML/text inventories, current-source name/signature
mapping and detailed checks of scientific/dispatch entry points. This is **not**
a re-derivation of every Numerical Recipes helper algorithm. No browser/JavaScript
renderer was available: live visual appearance and MathJax execution are not claimed
verified. However, HTML itself proves specific equation corruption: an unescaped
absolute-value pipe splits the k-normalization table and loses its d_ref expression;
raw dollar-math underscores become `<em ...>` elements and strip TeX syntax.
Loading MathJax cannot restore syntax already removed by Markdown processing.

Raw downloaded HTML, probe outputs, logs, package environment and rendered docs
remain in a temporary audit directory. Persistent URL/status/link/hash inventories
are in [PLUME_DOCUMENTATION_SITE_INVENTORY.json](PLUME_DOCUMENTATION_SITE_INVENTORY.json).
The main [JSON](PLUME_DOCUMENTATION_OUTPUT_AUDIT.json) stores parameter definitions,
capabilities, schemas, reader checks, all26example outcomes, issues and roadmap.
External third-party links were not exhaustively fetched. In particular, an
access-restricted scientific publisher/ADS URL must not be called dead merely
because a request would require credentials.

## 3. Current capability matrix

Legend: **SUPPORTED** means the documented core path is exercised; it does not
promise every supplementary diagnostic, branch, collision operator or conditioning
regime. **PARTIALLY SUPPORTED** requires the listed restrictions. **UNSUPPORTED**
means forbidden or the evaluated path fails, even if input parsing/process exit
accepts it. **NOT VALIDATED** is used for scientific accuracy not established by
the available checks, especially remaining finite-rate off-diagonal/FPC questions.

'''
    profiles=['kinetic_rho_ref','kinetic_d_ref','fully_cold_d_ref','mixed_warm_reference_d_ref','mixed_cold_reference_d_ref']
    modes=list(dict.fromkeys(x['mode'] for x in d['capability_matrix']))
    text+=table(['Workflow','Kinetic rho','Kinetic d','Cold d','Warm-ref mixed d','Cold-ref mixed d'],[[mode]+[next(x['status'] for x in d['capability_matrix'] if x['mode']==mode and x['plasma']==pr) for pr in profiles] for mode in modes])
    text+='''
All cold-containing plasmas require collisionless d_ref. All wavevector scans
exclude kpar=0 at written **and hidden** points; fully cold perpendicular
single-point O/X dispersion is separate. Cold/mixed exact-parallel branches use
the analytic dispatcher; all-warm exact-parallel currently does not. Nonzero
oblique warm calculations are the verified general thermal path. Charge neutrality,
current balance and reference drift checks remain warnings; no automatic density,
velocity or frame adjustment occurs.

Diagnostics must be listed separately:

| Diagnostic | Kinetic-only | Cold-only | Mixed, either reference |
|---|---|---|---|
| E/B, density, velocity | Available for valid Ex=1 mode; all-warm geometry caveats above | Valid fluid response, subject to gauge/poles | Full-plasma fields; per-species kinetic/fluid moments |
| Species total + seven selected heating channels | Corrected collisionless or explicit historical legacy | Kinetic heating unavailable/NaN | Corrected warm species only; cold slots NaN |
| Formal complete wave energy | Corrected real-frequency weak-mode diagnostic | Includes cold species derivatives | Includes **every** physical species derivative |
| Finite-collision-rate corrected heating | Unsupported | Unsupported | Unsupported |
| FPC | Existing gyro/Cartesian thermal paths, with reader/accuracy caveats | Unsupported | Unsupported, including warm members of mixed plasma |

2D eigen/heating use **scan_input_2** flags; susceptibility tensor output
requires tensor=true in **both** scan inputs. Map eigen/tensor flags are
accepted but ignored; corrected map heating is rejected, historical map heating
ignored. Plasma-parameter/radial corrected heating is rejected; legacy kinetic
support does not imply modern mixed support. Radial option5 remains experimental.

## 4. Input inventory and validation

There are **89 distinct case-insensitive user-facing names in14namelist definitions**:
params, species, maps, guess, scan_input, fpc, radial_input and seven alternative
k_range layouts. The JSON has every name's spelling, type, explicit/default-or-
required status, allowed values, definition/units, applicability, dependencies,
validation caveats, source line, documentation matches and example/test fixture.
A mention-match is not proof of adequate documentation.

'''
    text+=table(['Namelist','Fields','Initialization / requirements'],[[g['group'],', '.join(g['parameters']),'Indexed species/guess/scan groups required as configured; see JSON.'] for g in d['namelist_definitions']])
    text+='''
Safe defaults explicitly present in code include rho_ref, faddeeva, adaptive,
bessel_tol1e-10, nmax512, nfixed128, corrected heating, nscan0, collision_type0,
Kn1e30, low_ntrue, mapnr/ni128, positive_rootsfalse, radialsubstep1, and FPC
electric selector0. Scan eigen/heating/tensor reset false per group. Many central
values **have no assigned safe default**: beta/k/velocity/counts/root guesses,
map bounds/logflags, species parameters and FPC grid bounds/spacing must be
supplied. Do not turn a compiler's zero-initialized storage into a public default.

Validation is uneven. Modern selector enums and Bessel ranges are explicit;
cold combinations/geometries are guarded. Counts, finite/global/warm density,
Kn/grid parameters and reference self-ratios lack comprehensive early checks.
The species1 mu warning is one-sided/nonterminal; physical self-ratios should
be explained rather than claiming enforcement that does not exist. This is
separate from warning-only equilibrium neutrality/current behavior.

Internal variables are not inputs merely because they are public in vars:
solver_*, nu_ns, auxiliary b_star, nbesmax15, root tolerance/iteration limits,
computemoment and Sokhotski-Plemelj constants are source/internal controls,
not namelist parameters. Deprecated new_low_n is no longer accepted; heat is
not an alias for the input heating, and va_ref_c/b_star are internal API concepts.
In particular, adaptive Bessels for susceptibility
do not make FPC distribution harmonics adaptive. No controls were renamed here. Reference self-ratios and the distinction between
a reference population and aggregate species density need an explicit semantic
decision before adding new hard validation; this audit does not mandate D1=1
or introduce neutrality/current constraints.

## 5. Normalization and scientific diagnostic definitions

Let physical reference R be species1, with alpha_R=Tperp,R/Tparallel,R,
mu_s=mR/ms, Q_s=qR/qs signed, D_s=ns/nR, u_s=Us/vA,R. Then
rho_R=v_tperp,R/abs(Omega_R), d_R=vA,R/Omega_R, and normalized w=omega/Omega_R.
The ratio includes **reference anisotropy**:

    k_rho = k_d * sqrt(beta_R*alpha_R)
    vtp_rho = sqrt(beta_R)*a, a=vA,R/c=vtp_d
    tau_rho,s = beta_R/tau_d,s
    tau_d,s = 2*Tparallel,s/(mR*vA,R^2)
    v_tparallel,s/vA,R = sqrt(mu_s*tau_d,s)
    beta_parallel,s = D_s*tau_d,s
    Omega_s/Omega_R = mu_s/Q_s
    omega_ps^2/Omega_R^2 = D_s*mu_s/(Q_s^2*a^2)
    delta_s = w - kpar_d*u_s

Drifts are Alfven-normalized in both conventions, not thermal-normalized. Cold
species use exactly tau_d=0, ignore anisotropy placeholders and require no
thermal scale. Cold reference requires beta_R=0 and can coexist with warm
species. The internal auxiliary scale is selected from the first warm tau_d;
alpha*=1 is a **normalization** anisotropy, not a replacement physical alpha_s.
It adds no particle population, pressure or energy. Values like beta*,V*,K*,t*
belong in a developer explanation, not the physical input table.

Physical equivalence is demonstrated by the existing normalization_rho_ref.in /
normalization_d_ref.in pair (beta.4,alpha_R1.5,a1e-4,kd_perp.05,kd_par.06).
For scans, convert each point: a beta-only d_ref scan holds native a,kd and
other tau_d fixed, whereas a beta-only rho scan holds legacy thermal/k/ratios
fixed. Same endpoints do not imply an equivalent interpolated/radial path.
Kn_d=Kn_rho*sqrt(beta_R*alpha_R) is required for equivalent collision frequency.

Fields retain E=deltaE/deltaEx, **Ex=1**, B=deltaB/deltaEx in Gaussian units.
Faraday gives w*B=(kappa cross E)/a; the rho form uses
K/(vtp_rho*sqrt(alpha_R)). Species velocity is deltaU/(c*deltaEx/B0), density
is (delta n/n0)/(deltaEx/B0), both complex. To reconstruct current retain the
equilibrium-drift density current. Ex=0/degenerate modes are unavailable in
this gauge, not assigned a new normalization. These meanings apply to
cold-reference cases without fictitious thermal reference units.

Corrected heating uses Qabs_s=Im(E†chi_s(wr)E), Ps_s=Qabs_s/H. It is signed
absorption/emission approximately contributing to -gamma/wr for weak isolated
modes, not dimensional power or irreversible thermodynamic heating. Channels
are n0 yy,yz,zy,zz,fulln0,transverse+1,transverse-1. The first4close to5;
5--7 need not sum to species total; summing all7double-counts n0. Higher
harmonics and omitted +/-1longitudinal/cross contributions are intentionally
not in those selected columns. Cold kinetic heating is NaN, not zero.

Corrected formal energy is

    H = HE + HB + sum_s Hs
    HE = abs(E)^2
    HB = abs(kappa cross E/(a*wr))^2
    Hs = E† d[wr*chi_H,s(wr)]/dwr E
    W = abs(deltaEx)^2/(16*pi) * H

All physical species contribute Hs, including cold species. Derivatives fix
physical k/plasma/reference state; real laboratory wr is used, not a species
Doppler frequency. Signed H fractions can be negative or exceed1; f_E_abs=
HE/(abs(HE)+abs(HB)+sum(abs(Hs))) is a noncancelling contribution measure,
**not** a physical energy fraction. Unresolved H gives unavailable signed
fractions; formal negative H is not automatically an error. Strong/reactive
applicability flags qualify the weak-mode interpretation, not dispersion roots.
Prior independent audits/tests validate these corrected definitions in their
stated regimes, not exact conservation for arbitrary complex-frequency modes.

Legacy reproduction is three independent choices: heating_conventionlegacy,
z_backendlegacy, bessel_modefixed128 (plus source revision for full historical
physics replay). Z+fixed128 alone reproduces only numerical dispersion choices;
it does not restore old heating or revert the collisional zz correction.
FPC still uses historical thermal/debug paths. Its phase-speed header literally
stores w*sqrt(beta)/Kpar; with perpendicular rho and alpha_R!=1 this differs
from physical w/kpar_d by sqrt(alpha_R). Use exact header definitions, not
unqualified Alfven-speed labels. FPC header/reader normalization needs a separate scientific audit; the
velocity-integral audit remains distinct from any future mixed-FPC implementation.

The dispersion API stores ordinary dimensionless species chi, with epsilon=I+sum(chi)
and Lambda=epsilon+n*n^T-n^2*I, n=kappa/(a*w) at complex mode w. Its determinant
is raw; warm-reference direct stdout separately multiplies it by solver_vtp^6.
Historical comments about multiplying the wave tensor by w^2 do not describe
the implemented current expressions. Output/diagnostic documentation must not
confuse chi, epsilon, Lambda, raw determinant and displayed determinant scaling.
The normalization formulas above use PLUME's positive reference-frequency
convention; signed species frequencies are mu/Q. Broader physical reference
conventions beyond the tested input ratios should not be advertised without review.

## 6. Page-by-page website and source review

The final appendix lists **every probed HTML URL**, including unavailable
expected current pages. P0/P1entries give corrections/destinations; generic
Numerical Recipes/API entries remain useful as signatures but are not a current
scientific user guide. The JSON also inventories all repository root Markdown
sources. Existing detailed mathematics should be retained in topical scientific
pages; implementation chronology/reports should be visibly archived.

The source/published comparison is unambiguous: public input/output pages omit
new defaults and still contain old heating indices, while current repository
input/output sources include corrected controls and companion schema. Local
normalization.md remains a historical development report with contradictory
old/new stage restrictions. The landing page and Readme are separate home-like
pages with different text; titles such as Readme obscure purpose. The JET guide
and notebook links resolve404 because workflow staging does not publish them.
API navigation foregrounds generic helpers while current explicit-state cold,
energy, moment and cache modules are absent from the deployed API.

## 7. Output schemas and downstream compatibility

The main JSON supplies filename patterns, exact order/formula, producer,
normalization, species/harmonic order, NaNs/statuses, readers and limitations
for every output family. Core rectangular layouts are:

'''
    text+=table(['Output','Columns / structure','Primary caution'],[[x['id'],x['columns'],x.get('notes','')] for x in d['output_schemas']])
    text+='''
For low_n=true scans, zero-based heat_start=18+8S with eigen, otherwise6;
Ps(s)=heat_start+s; channel(s,c)=heat_start+S+7s+c; parameter suffix starts
heat_start+8S. Species ordering is original and +1precedes-1. Modes have native
k/beta/vtp prefix; susceptibility tensors instead have only4prefix columns,
then six independent tensor components in **component-major** order across
species, then species-major parameters. They do not contain frequency columns.

Map rows are5columns and (nr+1)*(ni+1)nodes, not nr*ni samples. Roots are6+6S;
the map writer currently writes all refined candidates before final filters,
not just nroot_max. Example_map_par has15root rows while continuing4branches.
There is also a hard gamma<1selection in map search, distinct from positive_roots,
which refers to real frequency. Roots/selection/continuation are not synonymous.

Corrected .energy has27+15S numeric fields plus quoted source_file and a versioned
header. Record/workflow/scan/outer/point/mode and root_valid identify mode rows;
point is the **internal** integration index, not row number with nres>1. Existing
mode precision/columns are unchanged. Tests/actualoutputs verify all4,831completed
associations, including failed roots; failed hidden substeps warn but have no
written record. An interrupted process can flush the companion before a mode
buffer, so completion/row-count checks remain necessary. Use unique run/output
names; generic folder globs can collect stale modes from an earlier run.

Recheck of the ten previously suspected output issues:

| Issue | Current classification | Evidence |
|---|---|---|
| Seven-channel Python indexing | **Already fixed** in current/debug explicit readers; public docs stale |32markerconfigurations pass; old aliases remain debug/legacy |
| Heating-only starting offset | **Already fixed** current reader | start6 marker totals correct |
| Species parameter offsets | **Already fixed** current supported low_n=true reader | all species/flags marker checks |
| Historical Ew | **Partially corrected intentionally** | correctedNaN; legacy undefined and not H |
| low_n=false formatting | **Partially corrected intentionally** | corrected orphan removed; legacy retained; main reader4IndexErrors |
| Corrected .energy companion | **Implemented and readable** | version/schema,quotedpath/NaN tests and actualoutputs |
| Companion/root association | **Correct for completed audited runs** |4,831records; interrupted output not atomic |
| Cold kinetic-heating NaN | **Correctly retained** active sweep/energy readers | markers andactualcold/mixed examples |
| Failed-root scan output | **Explicit NaN/status and recovery** | actualmixed scan9rows incl1failed; priorindependenttests |
| Map diagnostic requests | **Still silently omitted for eigen/tensor/legacyheat** | isolatedoption4probes; correctedheat rejected |

FPC is not a rectangular dispersion-mode format: gyro files have labeled
thermal/grid headers and2Dvelocity matrices; Cartesian files have projection/3D
blocks. Debug .eigen lacks a species-parameter suffix, while .mom orders density,
velocities and j.Edebugvalues. Do not feed them to a scan reader based on a
similar-looking prefix. Radial files append radius **last**, after any historical Ew; generic reading
can silently lose radius and its total-only heat block is not the seven-channel
scan schema. Runtime fmt in radial_scan is the detailed old-layout authority;
that experimental family is not declared scientifically validated for arbitrary
models merely because a synthetic interpolation regression passes.

## 8. Python, wrapper and plotting findings

Current load_plume_sweep with explicit nspec/eigen/heating reads supported
low_n=true files correctly; load_plume_sweep_debug returns **column indices**,
not values. Autodetection is mathematically ambiguous and its debug shortcut
still mentions older lengths; explicit flags are the recommended immediate
workflow. Generic2D reading flattens blank-separated blocks and does not retain
outer/inner indices. Read each mode file separately; use .energy context for joins.
load_plume_energy correctly preserves unavailable/failed rows and quoted paths.

The wrapper's writeOut argument is wrongly derived from use_map (2of4tested
combinations wrong). load_from_file drops low_n,loggridg,nr,ni,tensor and all
sixCartesian FPCvelocity bounds on roundtrip. make_scan does not expose tensor,
though manual dict entries can. set_maps does not expose all map fields. Several
execution helpers mutate option/use_map/scan state, invoke shell commands,
append logs and do not check exit codes, so old files can be mistaken for a new
successful run. Document that they are root-working-directory helpers with
side effects rather than a general safe process API.

loadeigen/loadmoms are two-species debugging interfaces. loadeigen fabricates
Pe=0 and derived products omit proper conjugation/equilibrium drift density
current; do not use them for general species energy closure. Two definitions
of load_plume_sweep_orig exist; the latter shadows the former. FPC old readers
replaceNaNwith0 and label specie02electron regardless population identity.
Cartesian parsing recursively shifts header offsets on broad exceptions rather
than validating a schema. These need reader fixes, not new plasma physics.

The contours reader has normalization-aware labels/formulas but falsely infers
heat from trailing parameter columns when no heat was requested. Its thermal
velocity-reference plotting at cold beta0 is not a cold-reference resonance/VDF
implementation. SWIFT parses E/B/reference velocity but calls coordinatesrho;
it can reconstruct a numerical phase if units are consistently reinterpreted,
yet provides no general d_ref/unit contract. The gnuplot examples hard-code
2species, column indices, filenames andrho labels. linfpcplot's two-species power
plot assumes species1ions/species2electrons, ignores signed-negative traces in
log presentation and labels an approximate diagnostic as gamma/omega. None is
a generic arbitrary-mixture diagnostic plot without explicit restrictions.

Recommended immediate reading workflow: save input+revision+normalization/backend
settings with each run, select the actual file family, use explicit scan switches,
retain NaNs, read .energy statuses, associate by mode/point/source, and label
native units. A later unified typed schema reader should replace heuristic
inference without changing existing numeric column positions.

The JET guide's option paragraphs mix gyro and Cartesian descriptions despite
its headings: the executable is unambiguously **6=gyro,7=Cartesian**. Literal
raw correlations in code are -Qs*v_i^2/2 times Re(E_i*conj(df1/dv_i)), with
species thermal velocity coordinates. f0hat is exp[-(vpar-hatU)^2-vperp^2/alpha].
calc_fs1 has prefactor -i*sqrt(mu*tau_rho/beta)*(exbar/Qs)*f0hat/vtp_solver;
exbar=1 for gyro but sqrt(beta)/vtp_solver for Cartesian. Cartesian outputs sum
omitted axes with delv; gyro outputs integrate phi with the vperp Jacobian and
a delv/4 central-cell substitution. These literal definitions are preserved in
JSON. Their conversion to dimensional power and a common cross-coordinate
amplitude are **not independently established here**; do not guess a conversion
or infer it from rho/d normalization-equivalence tests. That is a concrete
scientific-normalization question for the separate FPC audit, not implemented
FPC development in this stage.

## 9. Existing examples: execution and suitability

Inventory scope is every26file in inputs/example plus tracked public examples;
unrelated untracked research input/output/plot trees were neither altered nor
interpreted. Regression input fixtures are separately identified in JSON.
Each attempt used exact input text in its own temporary data directory; no
existing user output was overwritten. Timeout budget20s/process, three workers.
Wall times are characterization, not production performance benchmarks.

'''
    text+=table(['Input','Model / norm / option','Result','Warnings / introduction suitability'],[[Path(e['path']).name,f"{e['composition']} / {e['effective']['normalization']} / {e['effective'].get('option')}", 'TIMEOUT (not proved broken)' if e['execution'].get('timeout') else 'exit '+str(e['execution'].get('returncode')),f"root failures={e['execution'].get('root_failures')}; energy warnings={e['execution'].get('energy_warnings')}; {e['new_user_suitability']}"] for e in d['examples']])
    text+='''
Wide historical maps/scans are valuable advanced stress cases, but fail as an
unqualified clean first-calculation narrative: example_map_par has123failed
internal roots, high_beta_par46; many written heating slots are qualified or
unresolved. The mixed cold-reference scan has one failed written point and
later recovers. The large2D example produced12,256companion records before
interruption, including82derivative failures; this is **not** evidence that a
completed2D calculation cannot work. The installed tutorial asks for a map plot,
but current plotter/example/cycle.sh selects example_kpar-E instead and requires
gnuplot, pdflatex and appropriate data. Tool installation/rendered plots were
not exercised here; reader/data contracts were.

Branch names are not certified from plausible frequencies or filenames. The
small cold/eigen/heating examples have prior independent branch/field benchmarks;
new whistler/instability lessons below require an explicit verified input and
reference before promotion. Near crossings/poles/stiff damping, continuation
retains its existing limitations; do not demand identical forward/back branch
selection through a real mode crossing.

## 10. Proposed curated primary collection (12lessons)

Existing files are recommendations for retention, not changes made by this
audit. Empty existing-file entries mean **proposed** examples, not files claimed
to exist/run. One lesson can link a small pair or companion variant rather than
creating redundant inputs. Numerical/backend stress cases remain secondary.

'''
    text+=table(['Slot / topic','Existing candidates','Expected benchmark / acceptance'],[[x['slot']+' '+x['topic'],', '.join(x['existing_files']) or 'Proposed from validated test fixture; not added',x['benchmark']] for x in d['curated_examples']])
    text+='''
## 11. Documentation structure, navigation and build quality

Retain scientific detail and useful existing pages, but separate current user
contracts from development chronology. Suggested navigation:

1. Home / model overview (one page, no ALPS capability carryover).
2. Installation (core vs tests vs plots vs docs).
3. Quick Start (one supplied root, bounded scan, inspect output/status).
4. Input Reference (generated/reviewed defaults and conditional namelists).
5. Normalization (physical reference, paired inputs, auxiliary distinction).
6. Species / Plasma Models (kinetic, cold, mixtures, drifts, geometry matrix).
7. Calculation Modes and Scans (continuation/failure/mode selection).
8. Eigenfunctions / Species Moments (units, Exgauge, unavailable modes).
9. Heating and Wave Energy (Ps,7channels,HE/HB/Hs,fractions,weak applicability).
10. Output Reference / Python Reading (schemas, flags, NaNs and join keys).
11. Examples (12curated lessons plus advanced/stress index).
12. Numerical Methods / Limitations / Troubleshooting / Legacy Replay.
13. Developer API and **versioned historical report archive**.

Current doc.yml runs on main push/PR, stages rootMarkdown into docs, rewrites
some .md links with sed, copies selected developer reports/images, invokesFORD,
and deploys docs_out to gh-pages on push. It does not run numerical tests despite
INSTALL claiming automatic correctness checks. There is no dedicated example,
link, output-schema or defaults-check workflow. No-token builds omit citingpapers.
The static staging list omits the modern heating report and JET/notebook/assets.
Relative .md links with anchors escape the suffix-only sed rewrite, and report-
local links to result JSON become broken deployed targets.

Reproducing that workflow in a **temporary clone** with unpinned latest
FORD7.0.13/requests2.31.0 and no ADSsecret succeeds, generating235HTMLpages,
23modules and149procedurepages, but leaves64missing target occurrences.
It warns that srd_dir,base_url and project-level ordered_subpage are ignored;
current FORD expects project_url and page-level ordered_subpage. Coloned rho_ref/
d_ref source comments are parsed as unknown metadata, hiding units; inactive
src/gk_disp_gf.f90 also warns during parsing. These are evidence that “buildexit0”
is insufficient. Actual current deployment uses no exposed toolchain lock/source
commit; no claim is made that this local latest-FORD build is exactly the live
pipeline's environment.

Recommended automated checks: pinned isolated docsbuild, unexpected warnings,
expected page/API manifest, internalpaths+anchors+assets, browser math smoke
screenshots, every public input parsed in a fresh data directory, bounded example
root/statusbenchmarks, versioned output schemas/unique markers, input-default
inventory diff, no-token bibliography fallback, and all numerical/legacy tests.
No check should silently replace failed data with0 or regenerate numerical goldens.

## 12. Prioritized findings and separately testable roadmap

Categories: **A** necessary code/reader corrections; **B** documentation;
**C** example improvements; **D** build/quality infrastructure. Scope,
dependencies, regression risk and independent implementation are stored per
issue in JSON. P0 means scientific misinterpretation/broken core workflow;
P1 means major missing capability/workflow; P2 organization/maintenance.

'''
    text+=table(['ID / category','Finding','Evidence and next action'],[[x['id']+' '+x['category'],x['title'],x['evidence']+' '+x['recommendation']] for x in issue])
    text+=table(['Stage','Scope','Acceptance / dependencies'],[[x['stage'],x['scope'],x['acceptance']+' '+x['depends']] for x in d['roadmap']])
    text+='''
These recommendations do not authorize a broad solver refactor. Map-origin and
all-warm parallel fixes need separate focused instructions/tests; FPC availability,
collisional tensor physics and perpendicular scans remain separate development
stages. Documentation can immediately describe observed restrictions without
waiting for every code correction. Preserve original golden numerical values,
file positions, species/harmonic order and explicit legacy controls.

## 13. Reproduction, limitations and audit-created files

All builds/examples ran in /tmp. The snapshot/copy setup and current git commit
are recorded in JSON. Reproduce by copying tracked files plus required test
measurement fixtures to an isolated checkout, linking/copying read-only Git
history for pinned regressions, then running make and the standard201test suite.
The audit scripts accept explicit destinations and never patch sources:

```sh
python3 tests/documentation_output_audit.py inventory --output /tmp/static.json
python3 tests/documentation_output_audit.py site --output /tmp/site.json --cache /tmp/site-html
MPLCONFIGDIR=/tmp/plume-mpl python3 tests/documentation_workflow_probes.py \
  --checkout /tmp/plume-copy --artifacts /tmp/workflows --output /tmp/workflows.json
MPLCONFIGDIR=/tmp/plume-mpl python3 tests/documentation_behavior_probes.py \
  --checkout /tmp/plume-copy --artifacts /tmp/behavior --output /tmp/behavior.json
MPLCONFIGDIR=/tmp/plume-mpl PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover \
  -s tests -p 'test_documentation_output_audit.py' -v
# Reproduce workflow build only in the copy; it rewrites/stages Markdown there.
# Optional: python3 tests/build_documentation_audit_report.py --artifacts <artifact-dir>
```

The added tests are **frozen characterization**, not desired future contracts
for confirmed bugs. Corrective tasks should update/retire those defect tests,
while retaining the audit evidence and original solver regression expectations.
There is no newly implemented physics or production correction in Stage5A.

Limitations: no JavaScript/browser rendering, exhaustive external-link verification,
production toolchain matrix or new independent branch identification for thousands
of stress-scan points; no inference that interrupted runs are complete. Finite-rate
Krook off-diagonal/FPC scientific questions remain NOT VALIDATED as previously
reported. Untouched research directories are excluded from a public-example
inventory. Published absent pages are documented as unavailable, not reviewed.

Only these audit artifacts were created:

'''
    text+='\n'.join('- `'+x+'`' for x in d['created_files'])+'\n\n## Appendix: every published/probed HTML page\n\n'
    text+=table(['URL','Existing description','Verified behavior','Correction / priority / destination'],[[x['url'],x['existing_description'],x['verified_behavior'],x['required_correction']+' '+x['severity']+' → '+x['recommended_destination']] for x in d['documentation_pages']])
    for old,new in {'14namelist':'14 namelist','32marker':'32 marker','4,831completed':'4,831 completed','all26example':'all 26 example','26file':'26 file','20s/process':'20 s/process','2of4':'2 of 4','4IndexErrors':'4 IndexErrors','all7double':'all seven double','32 unique-marker':'32 unique-marker','beforepromotion':'before promotion','231':'231','case-insensitive user-facing names in14':'case-insensitive user-facing names in 14','positive_rootsfalse':'positive_roots=false','collision_type0':'collision_type=0','low_ntrue':'low_n=true','Nr128':'Nr=128','readers4':'readers 4','32markerconfigurations':'32 marker configurations','7stride':'7 stride','14namelistdefinitions':'14 namelist definitions','NoNaN':'No NaN','23modules':'23 modules','149procedurepages':'149 procedure pages','235HTMLpages':'235 HTML pages','201test':'201 test','82derivative':'82 derivative','123failed':'123 failed','46root':'46 root','onefailed':'one failed','mixedcold':'mixed cold','25of26':'25 of 26','primary12':'primary 12','Noheader':'No header','FPCvelocity':'FPC velocity','scan_input_2controls':'scan_input_2 controls'}.items():text=text.replace(old,new)
    path.write_text(text)

if __name__=='__main__':main()
