[minxval,minind]=min(expval_tgt_pop(2,:));

id_pop(minind)
%id: 732


[xval_s, ind_s] = sort(expval_tgt_pop(2,:));
%id: 732   289    33    19    36    37    56   858   995   931
%original id:  1064 889 793 319 117 1067 1065 1068 49 1072

%  id_pop(ind_s(1:10))'
%     {'hugo/2022/07July/11/20'     } %negative vision kernel
%     {'hugo/2021/09September/21/27'} not found?
%     {'hugo/2021/03March/23/17'    } negative vision kernel
%     {'hugo/2021/03March/22/22'    } not found?
%     {'hugo/2021/03March/23/22'    } negative vision kernel despite positive tgt resp
%     {'hugo/2021/03March/23/26'    } negative vision kernel
%     {'hugo/2021/03March/26/26'    } negative vision kernel despite flat tgt resp
%     {'hugo/2022/08August/12/20'   } not found?
%     {'hugo/2023/01January/26/20'  } flat kernel despite negative tgt resp
%     {'hugo/2022/09September/12/2' } not found?

%% vision 
[xval_v, ind_v] = sort(expval_tgt_pop(2,:),'descend');
%original ID:    198   195   141   150   583   552   589   183   568   228
% {'hugo/2021/09September/07/31'}
% {'hugo/2021/09September/07/26'}
% {'hugo/2021/08August/24/27'   }
% {'hugo/2021/08August/25/27'   }
% {'hugo/2022/07July/18/15'     }
% {'hugo/2022/07July/08/1'      }
% {'hugo/2022/07July/19/15'     }
% {'hugo/2021/09September/03/26'}
% {'hugo/2022/07July/12/2'      }
% {'hugo/2021/09September/21/25'}

%% eye speed
[xval_es, ind_es] = sort(expval_tgt_pop(3,:),'descend');
%original ID:  567   643   711   642   646   708   619   641   549   668


[ves, ind_ves] = sort(expval_tgt_pop(2,:)-expval_tgt_pop(3,:),'descend');


id_pop(ind_ves(2)) %inhibitory eye pos kernel
id_pop(ind_ves(3)) %not found
id_pop(ind_ves(4)) %not found

id_pop(ind_ves(end))
