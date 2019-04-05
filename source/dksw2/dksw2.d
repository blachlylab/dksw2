module dksw2.dksw2;

import dksw2.ksw2;
import dksw2.kalloc;

struct KSW2_profile
{
    byte[25] score_mat;
    ubyte[] target;
    byte gapo,gape,gapo2,gape2;
    // int8_t a = 2, b = 4, q = 4, e = 2, q2 = 13, e2 = 1;
	// int c, i, pair = 1, w = -1, flag = 0, rep = 1, zdrop = -1, no_kalloc = 0;
    int w=-1;
    int zdrop;
    int flag=0;
    this(string target,byte match=2,byte mismatch=4,byte gapo=4,byte gape=2,byte gapo2=13,byte gape2=1,int zdrop=-1){
        ksw_gen_simple_mat(5,score_mat.ptr,match,mismatch);
        this.target.length=target.length;
        for(auto i=0;i<target.length;i++){
            this.target[i]=seq_nt4_table[target[i]];
        }
        this.gapo=gapo;
        this.gapo2=gapo2;
        this.gape=gape;
        this.gape2=gape2;
        this.zdrop=zdrop;
    }
    KSW2_align_result align_query(string query){
        // void * km =km_init;
        ubyte[] q;
        for(auto i=0;i<query.length;i++){
            q[i]=seq_nt4_table[query[i]];
        }
        KSW2_align_result res;
        ksw_extd2_sse(null,cast(int)q.length,q.ptr,cast(int)target.length,target.ptr,5,score_mat.ptr,gapo,gape,gapo2,gape2,w,zdrop,0,flag,&res.ez);
        return res;
    }
}

struct KSW2_align_result
{
    ksw_extz_t ez;
    
}