function spm_TOM
% TOM Toolbox wrapper to call TOM functions
%_______________________________________________________________________
% Christian Gaser
% $Id$

SPMid = spm('FnBanner',mfilename,'v1.04');
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Template-O-Matic');
spm_help('!ContextHelp',mfilename);
spm_help('!Disp','TOM.man','',Fgraph,'      Template-O-Matic toolbox for SPM8');

fig = spm_figure('GetWin','Interactive');
h0  = uimenu(fig,...
	'Label',	'TOM',...
	'Separator',	'on',...
	'Tag',		'TOM',...
	'HandleVisibility','on');
h1  = uimenu(h0,...
	'Label',	'Estimate regression parameters',...
	'Separator',	'off',...
	'Tag',		'Estimate regression parameters',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.tom.estimate'');',...
	'HandleVisibility','on');
h2  = uimenu(h0,...
	'Label',	'Create new template',...
	'Separator',	'off',...
	'Tag',		'Create new template',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.tom.create'');',...
	'HandleVisibility','on');
h3  = uimenu(h0,...
	'Label',	'Print information about TOM.mat',...
	'Separator',	'off',...
	'Tag',		'Print information about TOM.mat',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.tom.getinfo'');',...
	'HandleVisibility','on');
