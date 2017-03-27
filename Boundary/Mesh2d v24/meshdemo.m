function meshdemo

% Demo function for mesh2d.
%
% Feel free to "borrow" any of the geometries for your own use.
%
% Example:
%
%   meshdemo;       % Runs the demos
%
% Darren Engwirda - 2006

clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

answer = input(['This is a demo function for mesh2d. \n'                                                     ...
                '\n'                                                                                         ...
                'Several example meshes are shown, starting with some simple examples and progressing to \n' ...
                'the CFD-like applications for which the function was designed. \n'                          ...
                '\n'                                                                                         ...
                'The following is a simple mesh in a circle. Continue?? [y/n] \n'],'s');

if ~strcmp(answer,'y')
    return
end

% Geometry
dtheta = pi/12;
theta  = (-pi:dtheta:(pi-dtheta))';
node   = [cos(theta) sin(theta)];

% Make mesh
[p,t] = mesh2d(node);


answer = input(['The element size function is generated automatically to try to adequately resolve the geometry. \n'    ...
                '\n'                                                                                                    ...
                'This means that the mesh size is related to the length of the line segments used to define the \n'     ...
                'geometry. The following example is the same as the last, but with more lines used to represent the \n' ...
                'circle. Continue?? [y/n] \n'],'s');

if ~strcmp(answer,'y')
    return
end

% Geometry
dtheta = pi/75;
theta  = (-pi:dtheta:(pi-dtheta))';
node   = [cos(theta) sin(theta)];

% Make mesh
[p,t] = mesh2d(node);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

answer = input(['It is often necessary to specify the element size in some locations and the following example \n' ...
                'illustrates the use of a user specified sizing.'                                                  ...
                '\n'                                                                                               ...
                'Continue [y/n] \n'],'s');

if ~strcmp(answer,'y')
    return
end

close all

node = [
    0   0
    1   0
    1   1
    0   1
        ];

hdata.fun = @hfun1;

[p,t] = mesh2d(node,[],hdata);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Sliver regions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

answer = input(['Mesh2d can now deal with very fine "sliver" geometry features \n' ...
                '\n'                                                               ...
                'Continue [y/n] \n'],'s');

if ~strcmp(answer,'y')
    return
end

close all



        node  = [0 0
                 3 0
                 3 3
                 0 3
                 0.1 1
                 0.11 1
                 0.11 2
                 0.1 2
                 ];
        cnect = [1 2
                 2 3
                 3 4
                 4 1
                 5 6
                 6 7
                 7 8
                 8 5
                 ];

        [p,t] = mesh2d(node,cnect);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Cylinder in crossflow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

answer = input(['The following is a mesh used to simulate the flow past a cylinder. \n'                ...
                '\n'                                                                                   ...
                'This example also shows how user specified and automatic size functions are combined' ...
                '\n'                                                                                   ...
                'Continue [y/n] \n'],'s');

if ~strcmp(answer,'y')
    return
end

close all

theta = 0:pi/100:(2*pi-pi/100);
x     = cos(theta)/2;
y     = sin(theta)/2;

node = [ x'  y'
        -5  -10
        25  -10
        25   10
        -5   10 
       ];

n = size(node,1)-4;

cnect = [(1:n-1)' (2:n)'
          n        1
          n+1      n+2
          n+2      n+3
          n+3      n+4
          n+4      n+1
        ]; 

hdata = [];
hdata.fun = @hfun2;

[p,t] = mesh2d(node,cnect,hdata);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Airfoil + flap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

answer = input(['The following is a mesh used to simulate the flow past an airfoil/flap configuration. \n' ...
                '\n'                                                                                       ...
                'Continue [y/n] \n'],'s');

if ~strcmp(answer,'y')
    return
end

close all

wing = [
 1.00003  0.00126  
 0.99730  0.00170  
 0.98914  0.00302  
 0.97563  0.00518  
 0.95693  0.00812  
 0.93324  0.01176  
 0.90482  0.01602  
 0.87197  0.02079  
 0.83506  0.02597  
 0.79449  0.03145  
 0.75070  0.03712  
 0.70417  0.04285  
 0.65541  0.04854  
 0.60496  0.05405  
 0.55335  0.05924  
 0.50117  0.06397  
 0.44897  0.06811  
 0.39733  0.07150  
 0.34681  0.07402  
 0.29796  0.07554  
 0.25131  0.07597  
 0.20738  0.07524  
 0.16604  0.07320  
 0.12732  0.06915  
 0.09230  0.06265  
 0.06203  0.05382  
 0.03730  0.04324  
 0.01865  0.03176  
 0.00628  0.02030  
 0.00015  0.00956  
 0.00000  0.00000  
 0.00533 -0.00792  
 0.01557 -0.01401  
 0.03029 -0.01870  
 0.04915 -0.02248  
 0.07195 -0.02586  
 0.09868 -0.02922  
 0.12954 -0.03282  
 0.16483 -0.03660  
 0.20483 -0.04016  
 0.24869 -0.04283  
 0.29531 -0.04446  
 0.34418 -0.04510  
 0.39476 -0.04482  
 0.44650 -0.04371  
 0.49883 -0.04188  
 0.55117 -0.03945  
 0.60296 -0.03655  
 0.65360 -0.03327  
 0.70257 -0.02975  
 0.74930 -0.02607  
 0.79330 -0.02235  
 0.83407 -0.01866  
 0.87118 -0.01512  
 0.90420 -0.01180  
 0.93279 -0.00880  
 0.95661 -0.00621  
 0.97543 -0.00410  
 0.98901 -0.00254  
 0.99722 -0.00158  
 0.99997 -0.00126   
 ];

flap = rotate(0.4*wing,10);
flap = move(flap,0.95,-0.1);
wing = rotate(wing,5);
wing = move(wing,0,0.05);

wall = [
      -1      -3
       4      -3
       4       3
      -1       3 
       ];

nwing = size(wing,1);
nflap = size(flap,1);
nwall = size(wall,1);

cwing = [(1:nwing-1)', (2:nwing)'; nwing, 1];
cflap = [(1:nflap-1)', (2:nflap)'; nflap, 1];
cwall = [(1:nwall-1)', (2:nwall)'; nwall, 1];

cnect = [
        cwing
        cflap+nwing
        cwall+nflap+nwing
        ];

node = [wing; flap; wall];

hdata = [];
hdata.edgeh = [(1:size(cnect,1)-4)', 0.005*ones(size(cnect,1)-4,1)];

options.dhmax = 0.25;

[p,t] = mesh2d(node,cnect,hdata,options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Lake Superior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

answer = input(['The following is a mesh of Lake Superior and is a standard test of mesh algorithms. \n' ...
                '\n'                                                                                     ...
                'This example uses the automatic size function only'                                     ...
                '\n'                                                                                     ...
                'Continue [y/n] \n'],'s');

if ~strcmp(answer,'y')
    return
end

close all

%       From:
%       Computer Solutions Europe AB
%       Bjornnasvagen 21
%       S-113 47 Stockholm
%       SWEDEN
%       Tel:    +46 8 15 30 22
%       Fax:    +46 8 15 76 35
%       E-mail: info@comsol.se
%       URL:    http://www.comsol.se/

p1=[
-8.9154147   1.6615920	
-8.8847064   1.5538653
-8.7270571   1.4396306
-8.6780329   1.3310253
-8.2763869   1.3642085
-7.7589437   1.5521208
-7.3399811   1.5866412
-7.2064877   1.7400111
-6.9871610   1.7838988
-6.9190874   1.6751018
-6.6594362   1.8237880
-6.4705078   2.0285109
-6.3982459   2.0257761
-6.1776033   2.1236043
-5.9683636   1.9041636
-6.0500631   1.6420917
-6.2751037   1.4382319
-6.2083449   1.2768143
-6.3593873   1.1233594
-6.4379782   0.9672913
-6.4389687   0.9408290
-6.3308569   0.9103191
-5.9216680   1.1606956
-5.8853530   1.1594465
-5.5306989   0.8827500
-5.4571118   0.9068911
-5.2776995   0.8218321
-4.6882658   1.0696109
-4.4325551   1.1157840
-3.9867787   1.5284986
-3.6232783   1.5733560
-3.5156848   1.5181578
-3.1151141   1.6162773
-2.6031886   1.9253444
-2.4201615   2.0814581
-2.0576784   2.1825510
-1.9469493   2.3929843
-1.6564457   2.6016249
-1.4385855   2.8113568
-1.1857180   2.9681802
-0.8965760   3.2839699
-0.4656316   3.4939341
-0.0715676   3.6520868
 0.3221563   3.5994536
 0.6089023   3.4943945
 0.6092879   3.3885375
 0.3942456   3.3878929
 0.2866788   3.4141384
 0.0717038   3.3345153
 0.0896864   3.2286668
 0.1076576   3.1757488
-0.1436793   3.0169900
-0.3235837   2.8584538
-0.4676923   2.7529346
-0.9391942   2.0668111
-1.0495193   1.7498742
-1.0882363   1.3531376
-1.0158452   1.3262275
-0.9424102   1.4845979
-0.6152380   1.7477533
-0.5061946   1.9061816
-0.2890742   2.0115316
-0.5066670   1.7473959
-0.6160022   1.5360392
-0.3979720   1.8000363
-0.1808403   1.8525961
 0.2894343   1.7998173
 0.5431103   1.6681146
 0.9058832   1.5373304
 0.9432856   1.3258124
 1.0534256   1.1147328
 1.1641909   0.8507940
 1.3109203   0.6930545
 1.3117251   0.5871977
 1.4214690   0.5351345
 1.9317215   0.5400872
 2.2226054   0.5965222
 2.3340522   0.4391798
 2.4073552   0.4137156
 2.4449497   0.3348350
 2.5909063   0.3369578
 2.6249835   0.4962917
 2.7343488   0.4979874
 2.8454551   0.3938967
 2.9366470   0.3954205
 2.9201944   0.2892562
 2.9931921   0.2905022
 3.3135431   0.7199108
 3.3863589   0.7213257
 3.9960329   1.1315466
 4.0705567   1.0538628
 4.6490632   1.1744760
 4.9039980   1.1550887
 5.4495532   1.1450567
 5.9850466   1.4278414
 6.1653805   1.4607402
 6.3476152   1.4409063
 6.6355959   1.5048353
 6.5309084   1.3948027
 6.4623697   1.2862027
 6.4803085   0.8098819
 6.4094979   0.7542354
 6.4100090   0.7012564
 6.4882671   0.5981826
 6.7068647   0.6065329
 6.8546905   0.5593290
 6.9296957   0.5093111
 7.3625916   0.6332439
 7.4195041   0.5826584
 7.4810256   0.4262426
 7.6622000   0.4606821
 7.7339418   0.4903887
 7.6927729   0.5946095
 7.6199031   0.5913803
 7.5788199   0.6956243
 7.4649955   0.7967085
 7.4604098   0.9025555
 7.6668095   1.1767640
 7.7709460   1.2875089
 7.7298813   1.3917132
 7.6707366   1.4951115
 7.6356734   1.4670289
 7.6246141   1.3074531
 7.5532118   1.2777792
 7.5280861   1.4357493
 7.5345391   1.7011603
 7.7131429   1.7621690
 7.8603274   1.7158544
 7.9228745   1.9308907
 7.7348899   2.0813402
 7.3713821   2.1182087
 7.2698318   1.9547751
 7.1885795   2.1633926
 7.1152657   2.1868122
 7.0441691   2.1573400
 6.9336956   2.2058041
 6.9676141   2.2602054
 7.0648074   2.5292990
 7.0897106   2.7954185
 7.2917362   3.1221593
 7.3230152   3.2295634
 7.2074219   3.4101256
 7.0987913   3.4319941
 7.0562590   3.5892427
 6.8349672   3.7390777
 6.5426013   3.8863685
 6.3538119   4.1440578
 6.6493349   4.7917061
 6.5385166   4.8932842
 6.7052993   5.1650313
 6.7009764   5.2708721
 5.7441174   5.2345817
 5.5686945   5.1755508
 4.7886070   5.1509481
 4.4636817   5.3536511
 4.2075761   5.6646354
 3.6075780   7.1858876
 3.3873644   7.6576363
 3.2804384   7.7613007
 3.1762847   7.7327191
 3.0381962   7.6506117
 2.8286881   7.6467634
 2.8598071   7.8590921
 2.6845531   7.9089847
 2.5468001   7.8008409
 2.4437671   7.6933492
 2.3390419   7.6917801
 2.1983122   7.7691861
 2.1616242   7.9010276
 2.0918985   7.9000964
 2.0235167   7.7933405
 1.9540682   7.7660064
 1.5705080   7.7353016
 1.5002116   7.7875578
 1.2547448   7.9442315
 0.7309650   8.1526824
 0.4522013   8.2575018
 0.1738657   8.3098811
-0.0694533   8.5215153
-0.1737495   8.4157382
-0.2084645   8.4422444
-0.3476150   8.3630957
-0.4518238   8.3898231
-0.5558115   8.4695476
-0.5215963   8.3106439
-0.5221188   8.1518584
-0.4879606   7.9400337
-0.2092654   7.8335660
-0.2097517   7.4630660
-0.5257651   7.0403597
-0.5604465   7.1463352
-0.7716285   6.9354925
-0.9473082   6.8835005
-0.9833619   6.7249256
-1.1593398   6.6731632
-1.2287946   6.7795396
-1.1574372   6.9378046
-0.6994015   7.4115277
-0.5933112   7.7286747
-0.6973173   7.8878837
-0.9759166   7.9422776
-1.1175557   7.6256239
-1.1194069   7.3609825
-1.2951661   7.2564438
-1.5483098   6.4117610
-1.9366363   6.3100700
-1.8637764   6.5209552
-1.7582839   6.5197737
-1.7202965   6.7840348
-1.5427314   6.9939696
-1.5417155   7.0998256
-1.5762336   7.1530937
-1.9985176   6.9989869
-2.3506579   6.8978706
-2.5645144   6.6893724
-2.5351899   6.3183276
-2.7562921   5.7395724
-3.7980607   5.0191439
-4.3373509   4.7676755
-4.9834948   4.5740477
-5.9256532   4.1282314
-7.0382285   3.1643631
-7.4517045   2.7580044
-8.2033903   2.1245258
-8.7235361   1.8639192
];

p2=[
-6.0222359   1.4026510
-5.9152018   1.5578835
-5.8726343   1.5299217
-5.7920161   1.5536639
-5.8345898   1.5816055
-5.8318773   1.6609928
-5.6635570   1.7083083
-5.5966615   1.7590827
-5.6808641   1.8413432
-5.7929717   1.7391469
-6.0818599   1.5637339
-6.0865686   1.4314222
];

p3=[
-5.7324455   1.9225485
-5.5683675   1.9594866
-5.5335991   1.9159672
-5.5318777   1.9688923
-5.4498160   2.0457126
-5.4914569   2.1000333
-5.7297704   2.0019357
-5.7677013   1.9502324
];

p4=[
-5.2706487   2.2254789
-5.2257103   2.2770640
-5.2196558   2.3563407
-5.1911650   2.4614151
-5.2924029   2.4539566
-5.3166868   2.3699543
-5.3193455   2.2852741
-5.3444930   2.1748120
];

p5=[
-2.6285090   4.8904087
-2.5231726   4.8093509
-2.4524965   4.7818143
-2.1664250   4.9101487
-1.8814086   4.9860466
-2.0588911   4.9881727
-2.0937130   5.0415487
-1.9150532   5.1452403
-1.7011793   5.2487523
-1.5239785   5.2470098
-1.1676824   5.5087407
-1.2022868   5.6148538
-0.9185046   5.7718016
-0.8828910   5.8245348
-0.9532116   5.8778612
-1.3075201   5.7215257
-1.5563954   5.5649116
-1.6276655   5.5126727
-1.7340945   5.4872998
-1.8405570   5.4619959
-2.0900008   5.3326514
-2.3756722   5.1777288
-2.5904890   5.0486276
];

p6=[
-0.1745905   7.6482740
-0.1048759   7.4629629
-0.0174793   7.4629295
-0.0349355   7.5687895
 0.0349007   7.7275752
 0.2267794   7.7806613
 0.3839078   7.7280333
 0.5583180   7.7550129
 0.6632222   7.7024853
 0.8201638   7.7296796
 0.7325480   7.8086476
 0.6273773   7.9405222
 0.4880418   7.9135694
 0.1393711   7.9922753
 0.0348428   7.9922181
-0.1743299   7.8864526
];

p7=[
 4.5964504   4.4039327
 4.8814214   4.4120893
 4.9866866   4.4681990
 5.1274771   4.5254771
 5.1614199   4.5795088
 5.0886198   4.6302255
 4.9818982   4.6269721
 4.8364965   4.7285905
 4.5535316   4.6675729
 4.3785074   4.5569131
 4.3813116   4.4510625
];


n1 = size(p1,1);
n2 = size(p2,1);
n3 = size(p3,1);
n4 = size(p4,1);
n5 = size(p5,1);
n6 = size(p6,1);
n7 = size(p7,1);

c1 = [(1:n1-1)', (2:n1)'; n1, 1];
c2 = [(1:n2-1)', (2:n2)'; n2, 1];
c3 = [(1:n3-1)', (2:n3)'; n3, 1];
c4 = [(1:n4-1)', (2:n4)'; n4, 1];
c5 = [(1:n5-1)', (2:n5)'; n5, 1];
c6 = [(1:n6-1)', (2:n6)'; n6, 1];
c7 = [(1:n7-1)', (2:n7)'; n7, 1];

node  = [p1; p2; p3; p4; p5; p6; p7];
cnect = [
        c1
        c2+n1
        c3+n2+n1
        c4+n3+n2+n1
        c5+n4+n3+n2+n1
        c6+n5+n4+n3+n2+n1
        c7+n6+n5+n4+n3+n2+n1
        ];

% Make mesh
[p,t] = mesh2d(node,cnect);


answer = input(['The following shows the influence of gradient limiting on the size function. \n' ...
                'The value dhmax is reduced to 0.1 \n'                                            ...
                '\n'                                                                              ...
                'Continue [y/n] \n'],'s');

if ~strcmp(answer,'y')
    return
end

close all

% Solver options
options.dhmax = 0.1;

% Make mesh
[p,t] = mesh2d(node,cnect,[],options);


% User defined size functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

answer = input(['The following example shows how the element size can be controlled using the  \n' ...
                'various settings in HDATA. \n'                                                    ...
                '\n'                                                                               ...
                'Continue [y/n] \n'],'s');

if ~strcmp(answer,'y')
    return
end

close all

hdata = [];
options.output = false;

node = [0 0; 1 0; 1 1; 0 1];              % Simple square example
% 
% [p,t] = mesh2d(node,[],[],options);       % Auto size fun only
% 
% figure
% subplot(2,2,1)
% patch('faces',t,'vertices',p,'facecolor','none','edgecolor','b');
% axis equal off;
% title('Automatic size fun only')
% 
% 
% hdata.hmax = 0.1;                         % Global size of 0.1
% 
% [p,t] = mesh2d(node,[],hdata,options);
% 
% subplot(2,2,3)
% patch('faces',t,'vertices',p,'facecolor','none','edgecolor','b');
% axis equal off;
% title('Global size fun with h0=0.1')


hdata.hmax = 0.1;                         % Global size of 0.1
hdata.edgeh = [1,0.05];                   % Boundary layer on bottom edge

[p,t] = mesh2d(node,[],hdata,options);

figure
subplot(2,2,2)
patch('faces',t,'vertices',p,'facecolor','none','edgecolor','b');
axis equal off;
title('Additional boundary layer function on bottom edge')

% hdata = [];
% hdata.fun = @hfun1;
% 
% [p,t] = mesh2d(node,[],hdata,options);
% 
% subplot(2,2,4)
% patch('faces',t,'vertices',p,'facecolor','none','edgecolor','b');
% axis equal off; hold on;
% title('User defined funcion centred at [0.25,0.75]')


% Refinement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

answer = input(['The following example shows how an existing mesh can be refined using the  \n' ...
                'REFINE function. This avoids doing expensive retriangulation. \n'              ...
                '\n'                                                                            ...
                'Continue [y/n] \n'],'s');

if ~strcmp(answer,'y')
    return
end

close all

hdata = [];
options.output = false;
hdata.hmax = 0.1;

node = [0 0; 1 0; 1 1; 0 1];              % Simple square example

[p,t] = mesh2d(node,[],hdata,options);    % Auto size fun only

figure
subplot(2,2,1)
patch('faces',t,'vertices',p,'facecolor','none','edgecolor','b');
axis equal off; hold on;
title('Uniform mesh from MESH2D')

node = [0.3,0.3; 0.7,0.3; 0.7,0.7; 0.3,0.7];

in = inpoly(p,node);
ti = sum(in(t),2)>0;
[p,t] = refine(p,t,ti);

subplot(2,2,2)
patch('faces',t,'vertices',p,'facecolor','none','edgecolor','b');
axis equal off; hold on;
title('Mesh refined in centre region using REFINE')

[p,t] = smoothmesh(p,t);

subplot(2,2,3)
patch('faces',t,'vertices',p,'facecolor','none','edgecolor','b');
axis equal off; hold on;
title('Smoothed mesh using SMOOTHMESH')

end      % meshdemo()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = move(p,xm,ym)

% Move a node set p by [xm,ym]

n = size(p,1);
p = p + [xm*ones(n,1), ym*ones(n,1)];

end      % move()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = rotate(p,A)

% Rotate a node set p by A degrees.

A = A*pi/180;
T = [ cos(A), sin(A)
     -sin(A), cos(A)];
p = (T*p')';

end      % rotate()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = hfun1(x,y)

% User defined size function for square

h = 0.01 + 0.1*sqrt( (x-0.25).^2+(y-0.75).^2 );

end      % hfun1()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = hfun2(x,y)

% User defined size function for cylinder

h1 = inf*ones(size(x));
in = (x>=0)&(x<=25)&(y>=-3)&(y<=3);
h1(in) = 0.2;

r = sqrt(x.^2+y.^2);
h2 = inf*h1;
h2(r<=3) = 0.02 + 0.05*r(r<=3);

h = min(h1,h2);

end      % hfun2()
