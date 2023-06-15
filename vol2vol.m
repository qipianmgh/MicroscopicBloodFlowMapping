function vol = vol2vol(mov,targ,R)
% vol = vol2vol(mov,targ,<R>)
%
% mov  = movable 3D vol
% targ = target 3D vol
%
% R is the transformation matrix
%
% Currently only uses nearest neighbor.
%
%


% This is modified from MRIvol2vol.m (See below) by David Boas 
%
% MRIvol2vol.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:09 $
%    $Revision: 1.3 $
%    $Modified by Qi Pian on Dec 2018
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

nwait = 6;
hwait = waitbar( 0, 'Transforming' );

vol = [];
if(nargin < 2 | nargin > 3)
  disp(sprintf('vol = vol2vol(mov,targ,<R>)\n'));
  return;
end

% Target vox to Mov vox Matrix
Vt2m = R;

nct = size(targ,1);
nrt = size(targ,2);
nst = size(targ,3);
nvt = nct*nrt*nst;


waitbar( 1/nwait, hwait, 'Transforming: computing transformed indices...');
fprintf('Computing transformed indices ... ');tic;
[tc tr ts] = meshgrid([1:nct],[1:nrt],[1:nst]);
tcrs = [tc(:) tr(:) ts(:) ones(nvt,1)]';
tcrs_trans = round(Vt2m*tcrs);
fprintf(' ... done %g\n',toc);

ncm = size(mov,1);
nrm = size(mov,2);
nsm = size(mov,3);
nvm = ncm*nrm*nsm;

waitbar( 2/nwait, hwait, 'Transforming: getting ok...');
fprintf('Getting ok ... ');tic
mc = tcrs_trans(1,:);
mr = tcrs_trans(2,:);
ms = tcrs_trans(3,:);
indok = find(mc >= 1 & mc <= ncm & ...
	     mr >= 1 & mr <= nrm & ...
	     ms >= 1 & ms <= nsm);
fprintf(' ... done %g\n',toc);
nok = length(indok);
fprintf('nok = %d\n',nok);

waitbar( 3/nwait, hwait, 'Transforming: getting tind of OCT...');
fprintf('Getting tind of OCT... ');tic
mc = mc(indok);
mr = mr(indok);
ms = ms(indok);
mind = sub2ind([ncm nrm nsm],mc,mr,ms);
fprintf(' ... done %g\n',toc);

waitbar( 4/nwait, hwait, 'Transforming: getting mind...');
fprintf('Getting mind ... ');tic
tc = tc(indok);
tr = tr(indok);
ts = ts(indok);
tind = sub2ind([nct nrt nst],tc,tr,ts);
fprintf(' ... done %g\n',toc);

waitbar( 5/nwait, hwait, 'Transforming: resampling...');
fprintf('Resampling ... ');tic
bg_level=min(mov(mind)); % set bg of transformed OCT image Qi added 08142018
vol = ones(nct,nrt,nst,1)*bg_level;
vol(tind) = mov(mind);
fprintf(' ... done %g\n',toc);

close(hwait)

return;



