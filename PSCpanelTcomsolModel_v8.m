function [model,TcellAve,TcellAveTop,TcellAveBack,TmodAve,TcellMax,Conv_AveTop,Conv_AveBack,Conv_AveSide,Rad_AveTop,Rad_AveBack,Rad_AveSide,CutLineCells,CutLineTop,CutLineBack,CutLineCellsX,CutLineTopX,CutLineBackX]=PSCpanelTcomsolModel_v8(Tamb,v_wind,cellDistance,Qcell,pathFolder,Tcell_init)

% Model exported on May 16 2024, 11:03 by COMSOL 6.2.0.290.

% min mesh size 0.02
% max mesh size 0.2

% Import comsol packages
import com.comsol.model.*
import com.comsol.model.util.*

% Create comsol model
model = ModelUtil.create('Model');

model.modelPath(pathFolder)

model.param.set('Tamb', [num2str(Tamb), '[degC]'], 'Ambient temperature');
model.param.set('v_wind', num2str(v_wind), 'Wind speed');
model.param.set('l', [num2str(cellDistance), '[cm]'], 'Distance between cells');
model.param.set('z_glass', '3000[um]', 'Thickness of glass');
model.param.set('z_eva', '500[um]', 'Thickness of EVA');
model.param.set('z_psc', '1.35[um]', 'Thickness of perovksite layer');
model.param.set('wp', 'z_glass+z_eva', 'Work plane z-coordinate');
model.param.set('x_cell', '15.6[cm]', 'Cell side length (x-direction)');
model.param.set('y_cell', '15.6[cm]', 'Cell side length (y-direction)');
model.param.set('ncx', '12', 'Number of cells x');
model.param.set('ncy', '6', 'Number of cells y');
model.param.set('x_panel', 'ncx*x_cell+(ncx+1)*l', 'Width (x-direction)');
model.param.set('y_panel', 'ncy*y_cell+(ncy+1)*l', 'Depth (y-direction)');
model.param.set('z_panel', '2*z_glass+z_eva', 'Height (z-direction)');
model.param.set('minm', '0.02[m]', 'Min mesh size');
model.param.set('maxm', '0.2[m]', 'Max mesh size');

% Comsol model component

model.component.create('comp1', true);

% Geometry
model.component('comp1').geom.create('geom1', 3);

model.result.table.create('tbl1', 'Table');
model.result.table.create('tbl2', 'Table');
model.result.table.create('tbl3', 'Table');
model.result.table.create('tbl4', 'Table');
model.result.table.create('tbl5', 'Table');
model.result.table.create('tbl6', 'Table');
model.result.table.create('tbl7', 'Table');
model.result.table.create('tbl8', 'Table');
model.result.table.create('tbl9', 'Table');
model.result.table.create('tbl10', 'Table');
model.result.table.create('tbl11', 'Table');
model.result.table.create('tbl12', 'Table');
model.result.table.create('tbl13', 'Table');
model.result.table.create('tbl14', 'Table');
model.result.table.create('tbl15', 'Table');

model.component('comp1').geom('geom1').create('blk1', 'Block');
model.component('comp1').geom('geom1').feature('blk1').label('Glass_bottom');
model.component('comp1').geom('geom1').feature('blk1').set('size', {'x_panel' 'y_panel' 'z_glass'});
model.component('comp1').geom('geom1').create('wp1', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp1').set('quickz', 'wp');
model.component('comp1').geom('geom1').feature('wp1').set('unite', true);
model.component('comp1').geom('geom1').feature('wp1').geom.create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r1').set('pos', {'l' 'l'});
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r1').set('size', {'x_cell' 'y_cell'});
model.component('comp1').geom('geom1').feature('wp1').geom.nodeGroup.create('grp1');
model.component('comp1').geom('geom1').feature('wp1').geom.nodeGroup('grp1').label('Cells');
model.component('comp1').geom('geom1').feature('wp1').geom.nodeGroup('grp1').placeAfter([]);
model.component('comp1').geom('geom1').feature('wp1').geom.nodeGroup('grp1').add('r1');
model.component('comp1').geom('geom1').create('arr1', 'Array');
model.component('comp1').geom('geom1').feature('arr1').label('Cell array');
model.component('comp1').geom('geom1').feature('arr1').set('fullsize', {'ncx' 'ncy' '1'});
model.component('comp1').geom('geom1').feature('arr1').set('displ', {'x_cell + l' 'y_cell + l' '0'});
model.component('comp1').geom('geom1').feature('arr1').selection('input').set({'wp1'});
model.component('comp1').geom('geom1').create('blk2', 'Block');
model.component('comp1').geom('geom1').feature('blk2').label('EVA');
model.component('comp1').geom('geom1').feature('blk2').set('pos', {'0' '0' 'z_glass'});
model.component('comp1').geom('geom1').feature('blk2').set('size', {'x_panel' 'y_panel' 'z_eva'});
model.component('comp1').geom('geom1').create('blk3', 'Block');
model.component('comp1').geom('geom1').feature('blk3').label('Glass_top');
model.component('comp1').geom('geom1').feature('blk3').set('pos', {'0' '0' 'z_glass + z_eva'});
model.component('comp1').geom('geom1').feature('blk3').set('size', {'x_panel' 'y_panel' 'z_glass'});
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

model.component('comp1').selection.create('sel1', 'Explicit');
model.component('comp1').selection('sel1').geom('geom1', 2);
model.component('comp1').selection('sel1').set([14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85]);
model.component('comp1').selection.create('sel2', 'Explicit');
model.component('comp1').selection('sel2').geom('geom1', 2);
model.component('comp1').selection('sel2').set([1 2 3 4 5 7 8 11 12 13 86 87 88]);
model.component('comp1').selection.create('sel3', 'Explicit');
model.component('comp1').selection('sel3').geom('geom1', 2);
model.component('comp1').selection('sel3').set([15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85]);
model.component('comp1').selection.create('sel4', 'Explicit');
model.component('comp1').selection('sel4').geom('geom1', 2);
model.component('comp1').selection('sel4').set([14]);
model.component('comp1').selection('sel1').label('Cells');
model.component('comp1').selection('sel2').label('Back surfaces');
model.component('comp1').selection('sel3').label('Cells-hotspot');
model.component('comp1').selection('sel4').label('Hotspot');

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material.create('mat2', 'Common');
model.component('comp1').material.create('mat3', 'Common');
model.component('comp1').material('mat1').selection.set([2]);
model.component('comp1').material('mat2').selection.set([1 3]);
model.component('comp1').material('mat3').selection.named('sel1');
model.component('comp1').material('mat3').propertyGroup.create('shell', 'Shell');

% Physics
model.component('comp1').physics.create('ht', 'HeatTransfer', 'geom1');
model.component('comp1').physics('ht').create('sls1', 'SolidLayeredShell', 2);
model.component('comp1').physics('ht').feature('sls1').selection.named('sel1');
model.component('comp1').physics('ht').feature('sls1').create('lhs1', 'LayeredHeatSource', 2);
model.component('comp1').physics('ht').feature('sls1').feature('lhs1').selection.named('sel1');
model.component('comp1').physics('ht').feature('sls1').create('lhs2', 'LayeredHeatSource', 2);
model.component('comp1').physics('ht').feature('sls1').feature('lhs2').selection.named('sel4');
model.component('comp1').physics('ht').create('hf1', 'HeatFluxBoundary', 2);
model.component('comp1').physics('ht').feature('hf1').selection.set([10]);
model.component('comp1').physics('ht').create('hf2', 'HeatFluxBoundary', 2);
model.component('comp1').physics('ht').feature('hf2').selection.named('sel2');
model.component('comp1').physics('ht').create('sar1', 'SurfaceToAmbientRadiation', 2);
model.component('comp1').physics('ht').feature('sar1').selection.set([10]);
model.component('comp1').physics('ht').create('sar2', 'SurfaceToAmbientRadiation', 2);
model.component('comp1').physics('ht').feature('sar2').selection.named('sel2');

model.component('comp1').physics('ht').feature('init1').set('Tinit', [num2str(Tcell_init),'[K]']);
model.component('comp1').physics('ht').feature('sls1').feature('iv1').set('Tinit_resistive', 'Tamb');
model.component('comp1').physics('ht').feature('sls1').feature('lhs1').set('Qs', [num2str(Qcell),'/z_psc']);
model.component('comp1').physics('ht').feature('sls1').feature('lhs2').set('Qs', 'if((x-x0)^2+(y-y0)^2<d^2/4,Q,0)');
model.component('comp1').physics('ht').feature('sls1').feature('lhs2').active(false);
model.component('comp1').physics('ht').feature('sls1').feature('lhs2').label('Hot spot');
model.component('comp1').physics('ht').feature('hf1').set('HeatFluxType', 'ConvectiveHeatFlux');
model.component('comp1').physics('ht').feature('hf1').set('h', '2.8+3*v_wind');
model.component('comp1').physics('ht').feature('hf1').set('Text', 'Tamb');
model.component('comp1').physics('ht').feature('hf1').label('Heat Flux (front surface)');
model.component('comp1').physics('ht').feature('hf2').set('HeatFluxType', 'ConvectiveHeatFlux');
model.component('comp1').physics('ht').feature('hf2').set('h', '(2.8+3*v_wind)/2');
model.component('comp1').physics('ht').feature('hf2').set('Text', 'Tamb');
model.component('comp1').physics('ht').feature('hf2').label('Heat Flux (back surface)');
model.component('comp1').physics('ht').feature('sar1').set('epsilon_rad_mat', 'userdef');
model.component('comp1').physics('ht').feature('sar1').set('epsilon_rad', 0.85);
model.component('comp1').physics('ht').feature('sar1').set('Tamb', 'Tamb');
model.component('comp1').physics('ht').feature('sar1').label('Surface-to-Ambient Radiation (top surface)');
model.component('comp1').physics('ht').feature('sar2').set('epsilon_rad_mat', 'userdef');
model.component('comp1').physics('ht').feature('sar2').set('epsilon_rad', 0.85);
model.component('comp1').physics('ht').feature('sar2').set('Tamb', 'Tamb');

% Mesh
model.component('comp1').mesh.create('mesh1');
model.component('comp1').mesh('mesh1').create('ftet1', 'FreeTet');
model.component('comp1').mesh('mesh1').create('ftet2', 'FreeTet');
model.component('comp1').mesh('mesh1').feature('ftet1').selection.geom('geom1', 3);
model.component('comp1').mesh('mesh1').feature('ftet1').selection.set([1 3]);
model.component('comp1').mesh('mesh1').feature('ftet1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftet2').selection.geom('geom1', 3);
model.component('comp1').mesh('mesh1').feature('ftet2').selection.set([2]);
model.component('comp1').mesh('mesh1').feature('ftet2').create('size1', 'Size');

model.result.table('tbl1').label('Surface Average Cells');
model.result.table('tbl1').comments('Surface Average Cells');
model.result.table('tbl2').label('Surface Average Module');
model.result.table('tbl2').comments('Volume Average 2');
model.result.table('tbl3').label('Surface Maximum Cells');
model.result.table('tbl3').comments('Surface Maximum 1');
model.result.table('tbl4').label('Surface Average, Top Surface');
model.result.table('tbl4').comments('Surface Average, Top Surface');
model.result.table('tbl5').label('Surface Average, Back Surface');
model.result.table('tbl5').comments('Surface Average, Back Surface');
model.result.table('tbl6').comments('Surface Average, Top Conv');
model.result.table('tbl7').comments('Surface Average, Back Conv');
model.result.table('tbl8').comments('Surface Average, Top Rad');
model.result.table('tbl9').comments('Surface Average, Back Rad');
model.result.table('tbl10').comments('Surface Average, Side Conv');
model.result.table('tbl11').comments('Surface Average, Side Conv');
model.result.table('tbl12').comments('Surface Average, Back Conv');
model.result.table('tbl13').comments('Surface Average, Side Conv');
model.result.table('tbl14').comments('Surface Average, Back Rad');
model.result.table('tbl15').comments('Surface Average, Side Rad');

% model.result.table('tbl6').comments('Surface Average 5');
% model.result.table('tbl7').comments('Surface Average 5');
% model.result.table('tbl8').comments('Surface Average 6');
% model.result.table('tbl9').comments('Surface Integration 1');
% model.result.table('tbl10').comments('Surface Average 7');
% model.result.table('tbl11').comments('Surface Average 8');

model.component('comp1').view('view1').set('transparency', true);
model.component('comp1').view('view2').axis.set('xmin', -0.5633813738822937);
model.component('comp1').view('view2').axis.set('xmax', 0.9398062825202942);
model.component('comp1').view('view2').axis.set('ymin', -0.30382171273231506);
model.component('comp1').view('view2').axis.set('ymax', 0.7960717082023621);

model.component('comp1').material('mat1').label('EVA');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', {'0.311' '0' '0' '0' '0.311' '0' '0' '0' '0.311'});
model.component('comp1').material('mat1').propertyGroup('def').set('density', '950');
model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', '2090');
model.component('comp1').material('mat2').label('Glass');
model.component('comp1').material('mat2').propertyGroup('def').set('thermalconductivity', {'2' '0' '0' '0' '2' '0' '0' '0' '2'});
model.component('comp1').material('mat2').propertyGroup('def').set('heatcapacity', '500');
model.component('comp1').material('mat2').propertyGroup('def').set('density', '2450');
model.component('comp1').material('mat3').propertyGroup('def').set('heatcapacity', '305.7');
model.component('comp1').material('mat3').propertyGroup('def').set('density', '5620.3');
model.component('comp1').material('mat3').propertyGroup('def').set('thermalconductivity', {'33' '0' '0' '0' '33' '0' '0' '0' '1'});
model.component('comp1').material('mat3').propertyGroup('shell').set('lth', 'z_psc');
model.component('comp1').material('mat3').propertyGroup('shell').set('lne', '5');

model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hmax', 'maxm');
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hmin', 'minm');
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hmax', 'maxm/2');
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hmin', 'minm/2');
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').run;

% Study
model.study.create('std1');
model.study('std1').create('stat', 'Stationary');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').create('i1', 'Iterative');
model.sol('sol1').feature('s1').create('d1', 'Direct');
model.sol('sol1').feature('s1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').create('so1', 'SOR');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').create('so1', 'SOR');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol1').feature('s1').feature.remove('fcDef');

% Results 
model.result.dataset.create('lshl1', 'LayeredMaterial');
model.result.dataset.create('cln1', 'CutLine3D');
model.result.dataset.create('cln2', 'CutLine3D');
model.result.dataset.create('cln3', 'CutLine3D');
model.result.dataset('lshl1').selection.geom('geom1', 2);
model.result.dataset('lshl1').selection.set([14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85]);
model.result.numerical.create('av1', 'AvSurface');
model.result.numerical.create('av2', 'AvVolume');
model.result.numerical.create('max1', 'MaxSurface');
model.result.numerical.create('av3', 'AvSurface');
model.result.numerical.create('av4', 'AvSurface');

model.result.numerical.create('int1', 'IntSurface');
model.result.numerical.create('int2', 'IntSurface');
model.result.numerical.create('int3', 'IntSurface');
model.result.numerical.create('int4', 'IntSurface');
model.result.numerical.create('int5', 'IntSurface');
model.result.numerical.create('int6', 'IntSurface');
model.result.numerical.create('int7', 'IntSurface');
model.result.numerical.create('int8', 'IntSurface');

model.result.numerical('av1').selection.named('sel1');
model.result.numerical('av2').selection.all;
model.result.numerical('max1').selection.named('sel1');
model.result.numerical('av3').selection.set([10]);
model.result.numerical('av4').selection.named('sel2');

model.result.numerical('int1').selection.set([10]);
model.result.numerical('int2').selection.named('sel2');
model.result.numerical('int3').selection.set([3]);
model.result.numerical('int4').selection.set([1 2 4 5 7 8 11 12 13 86 87 88]);
model.result.numerical('int5').selection.set([10]);
model.result.numerical('int6').selection.named('sel2');
model.result.numerical('int7').selection.set([3]);
model.result.numerical('int8').selection.set([1 2 4 5 7 8 11 12 13 86 87 88]);

% model.result.numerical('av5').selection.set([10]);
% model.result.numerical('av6').selection.named('sel2');
% model.result.numerical('av7').selection.set([10]);
% model.result.numerical('av8').selection.named('sel2');

model.result.create('pg1', 'PlotGroup3D');
model.result('pg1').selection.geom('geom1', 3);
model.result('pg1').selection.set([1 2 3]);
model.result('pg1').create('vol1', 'Volume');
model.result('pg1').create('vol2', 'Volume');
model.result('pg1').create('line1', 'Line');
model.result('pg1').feature('vol2').set('data', 'lshl1');
model.result('pg1').feature('line1').set('data', 'lshl1');
model.result('pg1').feature('line1').set('expr', '1');
model.result.export.create('tbl1', 'Table');
model.result.export.create('tbl2', 'Table');
model.result.export.create('tbl3', 'Table');
model.result.export.create('tbl4', 'Table');
model.result.export.create('tbl5', 'Table');
model.result.export.create('data1', 'Data');
model.result.export.create('data2', 'Data');
model.result.export.create('data3', 'Data');

model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').label('Compile Equations: Stationary');
model.sol('sol1').feature('v1').label('Dependent Variables 1.1');
model.sol('sol1').feature('s1').label('Stationary Solver 1.1');
model.sol('sol1').feature('s1').feature('dDef').label('Direct 2');
model.sol('sol1').feature('s1').feature('aDef').label('Advanced 1');
model.sol('sol1').feature('s1').feature('fc1').label('Fully Coupled 1.1');
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'i1');
model.sol('sol1').feature('s1').feature('fc1').set('initstep', 0.01);
model.sol('sol1').feature('s1').feature('fc1').set('minstep', 1.0E-6);
model.sol('sol1').feature('s1').feature('fc1').set('maxiter', 50);
model.sol('sol1').feature('s1').feature('fc1').set('termonres', false);
model.sol('sol1').feature('s1').feature('i1').label('AMG, heat transfer variables (ht)');
model.sol('sol1').feature('s1').feature('i1').set('nlinnormuse', true);
model.sol('sol1').feature('s1').feature('i1').feature('ilDef').label('Incomplete LU 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').label('Multigrid 1.1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('prefun', 'saamg');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('maxcoarsedof', 50000);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('saamgcompwise', true);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('usesmooth', false);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').label('Presmoother 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('soDef').label('SOR 2');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('so1').label('SOR 1.1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('so1').set('relax', 0.9);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').label('Postsmoother 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('soDef').label('SOR 2');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('so1').label('SOR 1.1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('so1').set('relax', 0.9);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').label('Coarse Solver 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('dDef').label('Direct 2');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('d1').label('Direct 1.1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol1').feature('s1').feature('d1').label('Direct, heat transfer variables (ht)');
model.sol('sol1').feature('s1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('s1').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol1').runAll;

model.result.dataset('cln1').label('Cut Line 3D Cells');
model.result.dataset('cln1').set('genpoints', {'0' '3*l+2.5*(y_cell)' 'wp'; 'x_panel' '3*l+2.5*(y_cell)' 'wp'});
model.result.dataset('cln2').label('Cut Line 3D Back Surface');
model.result.dataset('cln2').set('genpoints', {'0' '3*l+2.5*(y_cell)' '0'; 'x_panel' '3*l+2.5*(y_cell)' '0'});
model.result.dataset('cln3').label('Cut Line 3D Top Surface');
model.result.dataset('cln3').set('genpoints', {'0' '3*l+2.5*(y_cell)' '2*z_glass+z_eva'; 'x_panel' '3*l+2.5*(y_cell)' '2*z_glass+z_eva'});
model.result.numerical('av1').label('Surface Average Cells');
model.result.numerical('av1').set('table', 'tbl1');
model.result.numerical('av2').set('table', 'tbl2');
model.result.numerical('max1').set('table', 'tbl3');
model.result.numerical('av3').label('Surface Average, Top Surface');
model.result.numerical('av3').set('table', 'tbl4');
model.result.numerical('av4').label('Surface Average, Back Surface');
model.result.numerical('av4').set('table', 'tbl5');
% model.result.numerical('av5').label('Surface Average, Top Conv');
% model.result.numerical('av5').set('table', 'tbl7');
% model.result.numerical('av5').set('expr', {'ht.hf1.q0'});
% model.result.numerical('av5').set('unit', {'W/m^2'});
% model.result.numerical('av5').set('descr', {'Boundary convective heat flux'});
% model.result.numerical('av6').label('Surface Average, Back Conv');
% model.result.numerical('av6').set('table', 'tbl8');
% model.result.numerical('av6').set('expr', {'ht.hf2.q0'});
% model.result.numerical('av6').set('unit', {'W/m^2'});
% model.result.numerical('av6').set('descr', {'Boundary convective heat flux'});
% model.result.numerical('av7').label('Surface Average, Top Rad');
% model.result.numerical('av7').set('table', 'tbl10');
% model.result.numerical('av7').set('expr', {'ht.sar1.rflux'});
% model.result.numerical('av7').set('unit', {'W/m^2'});
% model.result.numerical('av7').set('descr', {'Radiative heat flux'});
% model.result.numerical('av8').label('Surface Average, Back Rad');
% model.result.numerical('av8').set('table', 'tbl11');
% model.result.numerical('av8').set('expr', {'ht.sar2.rflux'});
% model.result.numerical('av8').set('unit', {'W/m^2'});
% model.result.numerical('av8').set('descr', {'Radiative heat flux'});
% model.result.numerical('av1').setResult;
% model.result.numerical('av2').setResult;
% model.result.numerical('max1').setResult;
% model.result.numerical('av3').setResult;
% model.result.numerical('av4').setResult;
% model.result.numerical('av5').setResult;
% model.result.numerical('av6').setResult;
% model.result.numerical('av7').setResult;
% model.result.numerical('av8').setResult;
model.result.numerical('int1').label('Surface Int, Top Conv');
model.result.numerical('int1').set('table', 'tbl6');
model.result.numerical('int1').set('expr', {'ht.hf1.q0'});
model.result.numerical('int1').set('unit', {'W/m^2'});
model.result.numerical('int1').set('descr', {'Boundary convective heat flux'});
model.result.numerical('int2').label('Surface Int, BackAndSide Conv');
model.result.numerical('int2').set('table', 'tbl7');
model.result.numerical('int2').set('expr', {'ht.hf2.q0'});
model.result.numerical('int2').set('unit', {'W/m^2'});
model.result.numerical('int2').set('descr', {'Boundary convective heat flux'});
model.result.numerical('int3').label('Surface Average, Back Conv');
model.result.numerical('int3').set('table', 'tbl12');
model.result.numerical('int3').set('expr', {'ht.hf2.q0'});
model.result.numerical('int3').set('unit', {'W/m^2'});
model.result.numerical('int3').set('descr', {'Boundary convective heat flux'});
model.result.numerical('int4').label('Surface Average, Side Conv');
model.result.numerical('int4').set('table', 'tbl13');
model.result.numerical('int4').set('expr', {'ht.hf2.q0'});
model.result.numerical('int4').set('unit', {'W/m^2'});
model.result.numerical('int4').set('descr', {'Boundary convective heat flux'});
model.result.numerical('int5').label('Surface Average, Top Rad');
model.result.numerical('int5').set('table', 'tbl8');
model.result.numerical('int5').set('expr', {'ht.sar1.rflux'});
model.result.numerical('int5').set('unit', {'W/m^2'});
model.result.numerical('int5').set('descr', {'Radiative heat flux'});
model.result.numerical('int6').label('Surface Average, BackandSide Rad');
model.result.numerical('int6').set('table', 'tbl9');
model.result.numerical('int6').set('expr', {'ht.sar2.rflux'});
model.result.numerical('int6').set('unit', {'W/m^2'});
model.result.numerical('int6').set('descr', {'Radiative heat flux'});
model.result.numerical('int7').label('Surface Average, Back Rad');
model.result.numerical('int7').set('table', 'tbl14');
model.result.numerical('int7').set('expr', {'ht.sar2.rflux'});
model.result.numerical('int7').set('unit', {'W/m^2'});
model.result.numerical('int7').set('descr', {'Radiative heat flux'});
model.result.numerical('int8').label('Surface Average, Side Rad');
model.result.numerical('int8').set('table', 'tbl15');
model.result.numerical('int8').set('expr', {'ht.sar2.rflux'});
model.result.numerical('int8').set('unit', {'W/m^2'});
model.result.numerical('int8').set('descr', {'Radiative heat flux'});

model.result.numerical('av1').setResult;
model.result.numerical('av2').setResult;
model.result.numerical('max1').setResult;
model.result.numerical('av3').setResult;
model.result.numerical('av4').setResult;
model.result.numerical('int1').setResult;
model.result.numerical('int2').setResult;
model.result.numerical('int3').setResult;
model.result.numerical('int4').setResult;
model.result.numerical('int5').setResult;
model.result.numerical('int6').setResult;
model.result.numerical('int7').setResult;
model.result.numerical('int8').setResult;

% Extract results
TcellAve = model.result.numerical('av1').getReal();
TcellAveTop = model.result.numerical('av3').getReal();
TcellAveBack = model.result.numerical('av4').getReal();
TmodAve = model.result.numerical('av2').getReal();
TcellMax = model.result().numerical('max1').getReal();

% Conv_AveTop = model.result().numerical('av5').getReal();
% Conv_AveBack = model.result().numerical('av6').getReal();
% Rad_AveTop = model.result().numerical('av7').getReal();
% Rad_AveBack = model.result().numerical('av8').getReal();

% Convective Heat Fluxes
Conv_AveTop = model.result().numerical('int1').getReal();            %('int1').label('Surface Average, Top Conv');
Conv_AveBackAndSide = model.result().numerical('int2').getReal();    %('int2').label('Surface Average, BackAndSide Conv');
Conv_AveBack = model.result().numerical('int3').getReal();          %('int3').label('Surface Average, Back Conv');
Conv_AveSide = model.result().numerical('int4').getReal();          %('int4').label('Surface Average, Side Conv');

% Radiative Heat Fluxes 
Rad_AveTop = model.result().numerical('int5').getReal();         %('int5').label('Surface Average, Top Rad');
Rad_AveBackAndSide = model.result().numerical('int6').getReal(); %('int6').label('Surface Average, BackandSide Rad');
Rad_AveBack = model.result().numerical('int7').getReal();       %('int7').label('Surface Average, Back Rad');
Rad_AveSide = model.result().numerical('int8').getReal();       %('int8').label('Surface Average, Side Rad');

CutLineCells = mphinterp(model,'T','dataset','cln1');
CutLineTop = mphinterp(model,'T','dataset','cln3');
CutLineBack =  mphinterp(model,'T','dataset','cln2');

CutLineCellsX = mphinterp(model,'x','dataset','cln1');
CutLineTopX = mphinterp(model,'x','dataset','cln3');
CutLineBackX =  mphinterp(model,'x','dataset','cln2');

model.result('pg1').label('Temperature (ht)');
model.result('pg1').feature('vol1').label('Domain');
model.result('pg1').feature('vol1').set('unit', 'degC');
model.result('pg1').feature('vol1').set('colortable', 'HeatCameraLight');
model.result('pg1').feature('vol1').set('smooth', 'internal');
model.result('pg1').feature('vol1').set('resolution', 'normal');
model.result('pg1').feature('vol2').label('Layered Shell');
model.result('pg1').feature('vol2').set('unit', 'degC');
model.result('pg1').feature('vol2').set('titletype', 'none');
model.result('pg1').feature('vol2').set('smooth', 'internal');
model.result('pg1').feature('vol2').set('inheritplot', 'vol1');
model.result('pg1').feature('vol2').set('resolution', 'normal');
model.result('pg1').feature('line1').label('Layered Shell Edges');
model.result('pg1').feature('line1').set('titletype', 'none');
model.result('pg1').feature('line1').set('coloring', 'uniform');
model.result('pg1').feature('line1').set('color', 'fromtheme');
model.result('pg1').feature('line1').set('smooth', 'internal');
model.result('pg1').feature('line1').set('resolution', 'normal');

end
