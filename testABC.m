abc = DynSystem(@(t,x,e,eov)d_abc(t,x,e,eov), 3, [0,0]);
init = Grid(3, [1,2], [100,100], [0,2*pi; 0,2*pi], 1e-5);

ftl = FTLE(abc, init, [0,8], false);
imagesc(init.domain(:,1), init.domain(:,2), ftl);
colorbar;