function createArgs()
evalin('caller','syms px pz real;');
evalin('caller','py=sym("py","real");');
evalin('caller','syms vx vy vz real;');
evalin('caller','syms qw qx qy real;');
evalin('caller','qz=sym("qz","real");');
evalin('caller','syms bax bay baz real;');
evalin('caller','syms bgx bgy bgz real;');

evalin('caller','p=[px;py;pz];');
evalin('caller','v=[vx;vy;vz];');
evalin('caller','q=SymQuat(qw,qx,qy,qz);');
evalin('caller','ba=[bax;bay;baz];');
evalin('caller','bg=[bgx;bgy;bgz];');
evalin('caller','state=reshape([qw;qx;qy;qz;p;v;ba;bg],1,[]);');

evalin('caller','syms dx dy dz real;');
evalin('caller','delta=[dx; dy; dz];');
evalin('caller','delta_state=reshape([delta;p;v;ba;bg],1,[]);');


evalin('caller','g=[0;0;-9.81];');

%args=struct("dx", dx, "dy",dy,"dz",dz,"delta",delta,"delta_state",delta_state,"p",p,"v",v,"q",q,"ba",ba,"bg",bg,"state",state,"g",g,"px",px,"py",py,"pz",pz,"vx",vx,"vy",vy,"vz",vz,"qw",qw,"qx",qx,"qy",qy,"qz",qz,"bax",bax,"bay",bay,"baz",baz,"bgx",bgx,"bgy",bgy,"bgz",bgz);

%names = fieldnames(args);
%assignin('caller',"args",args);
%for i = 1:numel(names)
%    evalin('caller', join([names{i},'= args.',names{i},";"],""));
%end
end
