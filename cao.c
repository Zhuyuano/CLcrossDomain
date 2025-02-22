#include <pbc/pbc.h>
#include <pbc/pbc_test.h>



// 初始化元素
void init_elements(pairing_t pairing, element_t *P, element_t *BSN_A, element_t *s_A, element_t *P_pub_A, element_t *BSN_B, element_t *s_B, element_t *P_pub_B) {
    element_init_G1(*P, pairing);
	element_random(*P);
    element_init_G1(*BSN_A, pairing);
	element_init_G1(*P_pub_A, pairing);
    element_init_Zr(*s_A, pairing);
    element_random(*BSN_A);
    element_random(*s_A);
    element_pow_zn(*P_pub_A, *P, *s_A);
	
	element_init_G1(*BSN_B, pairing);
	element_init_G1(*P_pub_B, pairing);
    element_init_Zr(*s_B, pairing);
    element_random(*BSN_B);
    element_random(*s_B);
    element_pow_zn(*P_pub_B, *P, *s_B);
   
}
// 车辆注册函数
void vehicle_registration(pairing_t pairing, element_t P, element_t BSN, element_t P_pub, element_t *r, element_t *R, element_t *PID, element_t *SK, element_t *PK) {
    element_t ID,rPub,h0,h1,H1,H0;
    element_init_Zr(*r, pairing);
    element_init_G1(*R, pairing);
    element_init_Zr(ID, pairing);
	element_init_Zr(*PID, pairing);
	element_init_Zr(*SK, pairing);
	element_init_G1(*PK, pairing);
	
    element_random(*r);
    element_pow_zn(*R, P, *r);
	element_random(ID);
	
	element_init_G1(rPub, pairing);
    element_init_G1(H1, pairing);
    element_init_Zr(h1, pairing);
	element_init_G1(H0, pairing);
    element_init_Zr(h0, pairing);
	
	element_pow_zn(rPub, P_pub, *r);
	element_add(H1,BSN,rPub);
	element_from_hash(h1,"H1",1);
	element_pow_zn(H0,R,h1);
	element_from_hash(h0,"H0",0);
	element_add(*PID,ID,h0);

	element_random(*SK);
	element_pow_zn(*PK, P, *SK);
 
    element_clear(rPub);
	element_clear(h0);
	element_clear(H0);
	element_clear(h1);
	element_clear(H1);
	element_clear(ID);
}

// 车辆申请本地域证书
void vehicle_join(pairing_t pairing, element_t P, element_t PID, element_t R, element_t r, element_t s, element_t *t, element_t *c, element_t *Q) {
    element_t alpha,H2,f,H3,l,sc;
    element_init_Zr(alpha, pairing);
	element_init_Zr(H2, pairing);
	element_init_Zr(f, pairing);
	element_init_Zr(H3, pairing);
	element_init_Zr(l, pairing);
	element_init_Zr(sc, pairing);
	
    element_random(alpha);
	element_set1(H2);
	element_set1(H3);
	element_add(H2, H2, alpha);
	element_add(H2, H2, PID);
	element_from_hash(f,"H2",2);
    element_t beta, U;
    element_init_Zr(beta, pairing);
	element_init_G1(U, pairing);
	element_init_Zr(*t, pairing);
    element_init_G1(*Q, pairing);
	element_init_Zr(*c, pairing);


    element_random(beta);
    element_pow_zn(U, P, beta);
    element_add(*Q, U, R);
	element_add(H3, H3, PID);
	element_add(H3, H3, *Q);
	element_add(H3, H3, f);
	element_from_hash(*c,"H3",3);
    element_mul_zn(sc, s, *c);
    element_add(l, beta, sc);

    element_add(*t,r,l);
	
	element_clear(H2);
	element_clear(H3);
}

void authentication(pairing_t pairing, element_t P, element_t BSN_A, element_t P_pub_A, element_t BSN_B, element_t P_pub_B, element_t PID_i, element_t PID_j, element_t ti, element_t ci, element_t Qi, element_t tj, element_t cj, element_t Qj,element_t SK_j,element_t PK_j) {
    element_t n, y, N, E, Z, V, Y, Y1,nti,nci,temp,sigma,H4,H5,h5,H6,H7,Ki,Kj;
    element_init_Zr(n, pairing);
    element_init_Zr(y, pairing);
	element_init_G1(N, pairing);
	element_init_G1(E, pairing);
	element_init_G1(Z, pairing);
	element_init_G1(V, pairing);
	element_init_G1(Y, pairing);
	element_init_G1(Y1, pairing);
	element_init_Zr(nti, pairing);
	element_init_Zr(nci, pairing);
	element_init_Zr(temp, pairing);
	element_init_Zr(sigma, pairing);
	element_init_G1(H4, pairing);
	element_init_Zr(H5, pairing);
	element_init_Zr(h5, pairing);
	element_init_Zr(H6, pairing);
	element_init_Zr(H7, pairing);
	element_init_Zr(Ki, pairing);
	element_init_Zr(Kj, pairing);
	
    element_random(n);
    element_random(y);
	element_add(H4,BSN_A,P_pub_A);
	element_from_hash(N,"H4",4);
	
    element_mul_zn(nti, n, ti);
    element_pow_zn(E, P, nti);
    element_pow_zn(Z, Qi, n);
	element_mul_zn(nci, n, ci);
    element_pow_zn(V, P, nci);
	element_pow_zn(Y, P, y);
    element_pow_zn(Y1, PK_j, y); 


	element_add(H5,Y,Y1);
	element_add(H5,H5,PID_j);
	element_add(H5,H5,E);
	element_add(H5,H5,Z);
	element_add(H5,H5,V);
	element_from_hash(h5,"H5",5);
	element_mul_zn(temp, nti, h5);
    element_add(sigma, y, temp);


    // 验证 
	element_t temp1,temp2,temp3;
	element_init_GT(temp1, pairing);
	element_init_GT(temp2, pairing);
	element_init_GT(temp3, pairing);
    pairing_apply(temp1, E, P, pairing);
    pairing_apply(temp2, Z, P, pairing);
    pairing_apply(temp3, V, P_pub_A, pairing);


   
	// 使用 element_neg 来取负值
    element_t neg_SKj,Yj;
	element_init_G1(Yj, pairing);
	element_init_G1(neg_SKj, pairing);  // 初始化负数变量
	element_neg(neg_SKj, SK_j);
    element_pow_zn(Yj, Y1, neg_SKj);

    // 验证
	element_t left,right,temp5;
	element_init_G1(left, pairing);
	element_init_G1(right, pairing);
	element_init_G1(temp5, pairing);
	element_pow_zn(left,P,sigma);
	element_pow_zn(temp5,E,h5);
	element_add(right,Y,temp5);
	
	element_t w,W,temp6,temp7;
	element_init_Zr(w, pairing);
	element_init_G1(W, pairing);
	element_init_G1(temp6, pairing);
	element_init_G1(temp7, pairing);
	element_random(w);
	element_pow_zn(W,P,w);
	
	element_add(H6,PID_i,PID_j);
	element_add(H6,H6,Y);
	element_add(H6,H6,W);
	element_pow_zn(temp6,Y,w);
	element_add(H6,H6,temp6);
	element_from_hash(Kj,"H6",6);
	
	element_add(H7,PID_i,PID_j);
	element_add(H7,H7,Y);
	element_add(H7,H7,W);
	element_pow_zn(temp7,W,y);
	element_add(H7,H7,temp7);
	element_from_hash(Ki,"H7",7);
   
	element_clear(nti);
	element_clear(nci);
	element_clear(temp);
    element_clear(temp1);
	element_clear(temp2);
	element_clear(temp3);
	element_clear(left);
	element_clear(right);
	element_clear(temp5);
	element_clear(H4);
	element_clear(H5);
	element_clear(H6);
	element_clear(H7);
	element_clear(temp6);
	element_clear(temp7);
	element_clear(w);
	element_clear(W);
}

int main(int argc, char **argv) {
	double t0 = pbc_get_time();
    pairing_t pairing;
    element_t  PID_i, PID_j, BSN_A, BSN_B;
    element_t P, P_pub_A, P_pub_B, s_A, s_B, PK_i, PK_j, SK_i, SK_j, R_i, R_j,r_i,r_j, ti, ci, Qi, tj, cj, Qj;
	

    // 初始化配对
    pbc_demo_pairing_init(pairing, argc, argv);

    // 初始化元素
    init_elements(pairing, &P, &BSN_A, &s_A, &P_pub_A, &BSN_B, &s_B, &P_pub_B);
	
    // 初始化阶段
    double t1 = pbc_get_time();
    printf("Initialization Phase-----= %fs\n", t1 - t0);

    // 注册阶段
    double t2 = pbc_get_time();
    vehicle_registration(pairing, P, BSN_A, P_pub_A, &r_i, &R_i, &PID_i, &SK_i, &PK_i);
	vehicle_registration(pairing, P, BSN_B, P_pub_B, &r_j, &R_j, &PID_j, &SK_j, &PK_j);
    double t3 = pbc_get_time();
    printf("Registration Phase-----= %fs\n", t3 - t2);

    // 证书生成阶段 join
    double t4 = pbc_get_time();
    vehicle_join(pairing, P, PID_i, R_i, r_i, s_A, &ti, &ci, &Qi);
	vehicle_join(pairing, P, PID_j, R_j, r_j, s_B, &tj, &cj, &Qj);
    double t5 = pbc_get_time();
    printf("Join Phase-----= %fs\n", t5 - t4);
	
	// 身份认证阶段
    double t6 = pbc_get_time();
    authentication(pairing,P,BSN_A,P_pub_A, BSN_B, P_pub_B, PID_i, PID_j, ti, ci, Qi, tj, cj, Qj,SK_j,PK_j);
    double t7 = pbc_get_time();
    printf("Authentication Phase-----= %fs\n", t7 - t6);
	
    printf("Total time-----= %fs\n", t7 - t1);
    element_clear(P);
	element_clear(P_pub_A);
	element_clear(P_pub_B);
	element_clear(BSN_A);
	element_clear(BSN_B);
	element_clear(s_A);
	element_clear(s_B);
	pairing_clear(pairing);
	
}