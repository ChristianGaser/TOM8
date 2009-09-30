function varargout = cg_tom(cmd, job)
% Execution file for TOM
%_______________________________________________________________________
% Christian Gaser
% $Id$

linfun = inline('fprintf([''%-40s%s''],x,[repmat(sprintf(''\b''),1,40)])','x');

modi = {'GM','WM','CSF','BKG1','BKG2','BKG3','T1'};

if strcmpi(cmd,'estimate')
  % initialise output argument
  out = struct('tommat',{{''}},'betas',{{''}});
  GM   = job.data.GM;
  WM   = job.data.WM;
  CSF  = job.data.CSF;
  BKG1 = job.data.BKG1;
  BKG2 = job.data.BKG2;
  BKG3 = job.data.BKG3;
  T1   = job.data.T1;
  data = {GM, WM, CSF, BKG1, BKG2, BKG3, T1};

  age     = job.age.c;
  gender  = job.gender.c;
  agepoly = fieldnames(job.poly);
  odir = job.odir{1};
  covs = job.cov;
  n_covs  = length(covs);
  
  switch agepoly{1}
    case 'poly1', ageDegree = 1;
    case 'poly2', ageDegree = 2;
    case 'poly3', ageDegree = 3;
  end

  n_modi = zeros(1,7);
  for i=1:7
    n_modi(i) = length(data{i});
    % correct if entry is empty
    if strcmp(data{i},''), n_modi(i) = 0; end
  end

  % check for same file numbers
  if any(diff(n_modi(find(n_modi>0)))), error('Number of GM/WM/CSF/BKG?/T1 images differ.'); end
  % check for correct number of covariates
  if n_modi(1)~=size(age,1), error('Number of age values differ from number of images.'); end
  if n_modi(1)~=size(gender,1), error('Number of gender values differ from number of images.'); end
  for i = 1:n_covs
    if n_modi(1)~=size(covs(i).c,1)
      error(sprintf('Number of %s differ from number of images.',covs(i).cname));
    end
  end
  
  if (min(gender)<0) | (max(gender)>1)
    error('Valid values are 1 for males and 0 for females.');
  end

  % age regression with polynomial function
  X = zeros(n_modi(1),ageDegree+1);
  for j = 0:ageDegree, X(:,(j + 1)) = (age.^j) - X*(pinv(X)*(age.^j)); end

  % calculate correction for orthogonalization
  for j = 1:ageDegree, ageReg{j} = pinv(X(:,1:j))*(age.^j); end

  %-Ask about overwriting files from previous analyses...
  %-------------------------------------------------------------------
  if exist(fullfile(odir,'TOM.mat'),'file')
    str = {	'Current directory contains existing TOM.mat file:',...
            'Continuing will overwrite all files!'};
    if spm_input(str,1,'bd','stop|continue',[1,0],1,mfilename);
        fprintf('%-40s: %30s\n\n',...
            'Abort...   (existing TOM.mat file)',spm('time'));
        return
    end
  end

  % remove old files
  files = {'^T1_beta_.{2}\..{3}$','^GM_beta_.{2}\..{3}$','^WM_beta_.{2}\..{3}$',...
          '^CSF_beta_.{2}\..{3}$','^BKG.{1}_beta_.{2}\..{3}$','^TOM\..{3}$',...
          '^T1_Template_Age','^GM_Template_Age','^WM_Template_Age',...
          '^CSF_Template_Age','^BKG.{1}_Template_Age'};
 
  for i=1:length(files)
    j = spm_select('List',odir,files{i});
    for k=1:size(j,1)
      spm_unlink(fullfile(odir,deblank(j(k,:))));
    end
  end
  
  % extend design matrix with gender
  X = [X gender-mean(gender)];
  for i = 1:n_covs, X = [X covs(i).c-mean(covs(i).c)]; end

  n_beta = size(X,2);
  
  % save regression for orthogonalization
  out.tommat = {fullfile(odir,'TOM.mat')};
  save(out.tommat{1},'-V6', 'X', 'ageReg', 'age', 'gender', 'covs', 'n_modi');

  pKX = pinv(X);
  W = zeros(n_modi(1));
  for i = 1:n_modi(1), W(i,i) = 1; end
  W = sparse(W);

  TH = -Inf*ones(n_modi(1),1);

  % names of covariates
  for i = 1:ageDegree+1, betanames{i} = ['age' num2str(i-1)]; end
  betanames{i+1} = 'gender';
  for j = 1:n_covs, betanames{j+i+1} = covs(j).cname; end
    
  % go through all modi which are defined
  for i = find(n_modi>0)
  
	fprintf('Estimating %s...\n',modi{i});

    P = [];
    for j = 1:length(data{i})
      P = strvcat(P,data{i}{j});
    end
    V = spm_vol(P);

    % scale to a mean of 0.5 for t1 images
    if strcmp(modi{i},'T1')
      fprintf('Scale T1 images to mean of 0.5\n');
      for j = 1:n_modi(i)
        V(j).pinfo(1:2,:) = 0.5*V(j).pinfo(1:2,:)/spm_global(V(j));
      end
    end

    % estimate betas
    beta = calc_beta(V,pKX);

    dim = V(1).dim(1:3);

    VO = V(1);
    VO.dt = [spm_type('float32') spm_platform('bigend')];

    beta = reshape(beta,[dim n_beta]);
    
    % save beta images with prepending mod-name
    out.betas = cell(n_beta,1);
    for j=1:n_beta
      out.betas{j} = fullfile(odir,sprintf('%s_beta_%02d.img',modi{i},j-1));
      VO.fname = out.betas{j};
      VO.descrip = betanames{j};
      spm_write_vol(VO,beta(:,:,:,j));
    end

  end
  
elseif strcmpi(cmd,'create')
  % initialise output - does not pass on gif filename
  out = struct('templates',{{''}},'tpm',{{''}});
  
  newage = job.age.c;
  odir = job.odir{1};
  newcovs = job.cov;
  n_newcovs = length(newcovs);
  if isfield(job.newgender,'genderval')
    newgender = job.newgender.genderval;
  else
    newgender = [];
  end

  if isempty(newgender)
    is_gender = 0;
  else is_gender = 1; end
  
  % check for correct number of covariates
  if is_gender
    if (size(newgender,1)~=size(newage,1))
      error('Number of age values differ from number of gender values.');
    end
  end
  
  for i = 1:n_newcovs
      if size(newage,1)~=size(newcovs(i).c,1)
      error(sprintf('Number of %s differ from number of images.',newcovs(i).cname));
    end
  end
  
  load(job.tommat{:});
  idir = fileparts(job.tommat{:});

  ageDegree = size(ageReg,2);
   
  % check age range
  if max(newage) > max(age)
    error(sprintf('%3.2f is outside of the age range used for calculating the age regression (%3.2f..%3.2f years)\n',max(newage),min(age),max(age)));
  end
  if min(newage) < min(age)
    error(sprintf('%3.2f is outside of the age range used for calculating the age regression (%3.2f..%3.2f years)\n',min(newage),min(age),max(age)));
  end

  % check range of gender values
  if is_gender
    if (max(newgender) > 1) | (min(newgender) < 0)
      error('Gender should be defined only with 0 and 1.');
    end
  end
  
  % too many new covariates?
  if n_newcovs > length(covs)
    error('You have more covariates defined as used in the estimated model.');
  end

  % force a single mean value for the average approach
  if isfield(job.template,'average')
    newage = mean(newage);
    if is_gender, newgender = mean(newgender); end
    for i = 1:n_newcovs, newcovs(i).c = mean(newcovs(i).c); end
  end

  % vector of beta which will be used for template creation
  use_beta = 1:ageDegree + 2 + n_newcovs;
  
  % if no gender is defined do not use this beta
  if ~is_gender, use_beta(ageDegree+2) = []; end

  if isfield(job.gif,'gifslice')
    save_gif = 1;
    if size(newage,1) == 1
      disp('If you only have a single age value we cannot save animated gif-image.');
      save_gif = 0;
    end
  else save_gif = 0; end

  % go through all modi which are given
  for i = find(n_modi>0)
  
    % load beta images
    P = [];
    for j = use_beta
      P = strvcat(P,fullfile(idir,sprintf('%s_beta_%02d.img',modi{i},j-1)));
    end

    V = spm_vol(P);

    % calculate slice for gif-images
    if save_gif
      % voxelsize and origin
      vx =  sqrt(sum(V(1).mat(1:3,1:3).^2));
      Orig = V(1).mat\[0 0 0 1]';
      sl = job.gif.gifslice/vx(3)+Orig(3);
    end

    beta = spm_read_vols(V);

    VO = V(1);
    VO.dt = [spm_type('uint8') spm_platform('bigend')];
    VO.descrip = modi{i};
    VO.pinfo(1:3) = [1/255 0 0];

    template_sum{i} = zeros(size(beta(:,:,:,1)));
    n_newage = length(newage);
    
    if save_gif
      sz = size(beta(:,:,1,1));
      sz = [sz(2) sz(1) n_newage];
      img_array{i} = zeros(sz);
      gif_array = uint8([]);
    end
    
    for loop = 1:n_newage;  
      linfun(sprintf('Age: %3.2f years for %s',newage(loop), modi{i}));
      X = 1;
      for j = 1:ageDegree
        X = [X (newage(loop).^j) - X*ageReg{j}];
      end

      % extend design matrix with gender if defined
      if is_gender, X = [X newgender(loop)-mean(gender)]; end
      
      % extend covariates
      for j = 1:n_newcovs
        c = newcovs(j).c;
        X = [X c(loop)-mean(covs(j).c)];
      end

      % use offset
      template = beta(:,:,:,1);
      
      % and add weighted beta
      for j = 2:length(use_beta)
        template = template + beta(:,:,:,j)*X(j);
      end
      
      % prevent values outside the expected range 0..1
      template(find(template<0)) = 0;
      template(find(template>1)) = 1;

      template_sum{i} = template_sum{i} + template;
    
      % animated gif
      if save_gif
        gif_name = fullfile(odir,...
          [modi{i} '_Template_Age' num2str(newage(1)) '-' num2str(newage(end)) '_' num2str(job.gif.gifslice) 'mm.gif']);
        img = uint8(round(255*rot90(template(:,:,round(sl)))/max(template(:))));
        img_array{i}(:,:,loop) = img;
        sz = size(img,2);
        if n_newage > 1
          img(end-2:end,1:round(sz*(loop-1)/(n_newage-1))) = 255;
        end
        if loop==1
          imwrite(img,gif_name,'DelayTime',0.1,'LoopCount',Inf);
        else
          imwrite(img,gif_name,'WriteMode','Append','DelayTime',0.1);
        end
      end
            
    end
   
  end
  
  for i = find(n_modi>0)
  
    if save_gif
      gif_array = [gif_array img_array{i}];
    end
    
    % calculate mean of all ages
    if isfield(job.template,'matched')
      out.templates{i} = fullfile(odir,[modi{i} '_Template_Age' num2str(newage(1)) '-' num2str(newage(end)) '.img']);
    else
      out.templates{i} = fullfile(odir,[modi{i} '_Template_Age' num2str(newage(1)) '.img']);            
    end
    VO.fname = out.templates{i};
    template_sum{i} = template_sum{i}/n_newage;
    spm_write_vol(VO,template_sum{i});
    
    % also save mat-file to use theses images with older spm-versions
    % this idea is based on a script by Ged Ridgeway
    mat = VO.mat; M = mat;
    [pth fnm ext] = spm_fileparts(VO.fname);
    matfile = fullfile(pth, [fnm '.mat']);
    save(matfile, 'M', 'mat');    
  end
  
  % check that all 6 tissue classes were estimated to write TPM
  if ~any(n_modi(1:6)==0) 
    % calculate mean of all ages
    if isfield(job.template,'matched')
      out.tpm = fullfile(odir,['TPM_Age' num2str(newage(1)) '-' num2str(newage(end)) '.nii']);
    else
      out.tpm = fullfile(odir,['TPM_Age' num2str(newage(1)) '.nii']);            
    end
    VO.fname = out.tpm;

    TPM      = nifti;
    TPM.dat  = file_array(out.tpm,[size(template) 6], [spm_type('uint8') spm_platform('bigend')], 0, 1/255, 0);
    TPM.mat  = VO.mat;
    TPM.mat0 = VO.mat;
    create(TPM);
    
    all_sum = template_sum{1}+template_sum{2}+template_sum{3}+template_sum{4}+template_sum{5}+template_sum{6}+eps;
    for i=1:6
      template_sum{i} = template_sum{i}./all_sum;
      TPM.dat(:,:,:,1) = template_sum{i};
    end
  
  end

  % save combined gif image for all modalities
  if save_gif
    gif_name = fullfile(odir,['Template_Age' num2str(newage(1)) '-' num2str(newage(end)) ...
      '_' num2str(job.gif.gifslice) 'mm.gif']);
    sz = size(gif_array,2);
    for loop = 1:n_newage
      img = gif_array(:,:,loop);
      if n_newage > 1
        img(end-2:end,1:round(sz*(loop-1)/(n_newage-1))) = 255;
      end
      if loop==1
        imwrite(img,gif_name,'DelayTime',0.1,'LoopCount',Inf);
      else
        imwrite(img,gif_name,'WriteMode','Append','DelayTime',0.1);
      end
    end
  end
  
     
elseif strcmpi(cmd,'getinfo')
  % dummy output
  out = [];
  
  load(job.tommat{:});

  fprintf('\nInformation about %s\n\n',job.tommat{:});
  % print range of age and gender  
  fprintf('Age:\t%d values (%3.2f .. %3.2f years)\n',length(age),min(age),max(age));
  fprintf('Gender:\t%d females / %d males\n',sum(find(gender==0)>0),sum(find(gender==1)>0));
  
  % print information about covariates
  n_covs = length(covs);
  for i = 1:n_covs
    fprintf('Covariate %d:\t%s\n',i,covs(i).cname);
  end

  % print avaliable modi
  ind = find(n_modi>0);
  fprintf('Modi:\t');
  fprintf('%s ',modi{ind});
  fprintf('\n');
  
  % print order
  fprintf('Order of age regression: %d\n',length(ageReg));

end

% return output, if requested
if nargout > 0
    varargout{1} = out;
end
return

% -----------------------------------------------------------------------------------------------------

function beta = calc_beta(VY,pKX);

n = size(VY,1);
Y = zeros([prod(VY(1).dim(1:2)) n]);
beta = zeros([VY(1).dim(1:3) size(pKX,1)]);

spm_progress_bar('Init',VY(1).dim(3),'Computing','planes completed')
for j=1:VY(1).dim(3),

  M  = spm_matrix([0 0 j]);

  % Load slice j from all images
  for i=1:n
    tmp = spm_slice_vol(VY(i),M,VY(1).dim(1:2),[1 0]);
    Y(:,i) = tmp(:);
  end
  
  pXY = pKX*Y';
  pXY = reshape(pXY',[VY(1).dim(1:2) size(pKX,1)]);
  beta(:,:,j,:) = pXY;
  spm_progress_bar('Set',j);
end

spm_progress_bar('Clear');

return
