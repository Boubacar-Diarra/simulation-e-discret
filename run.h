
#ifndef PS_RUN_H
#define PS_RUN_H

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <cmath>
#include "mt.h"

using namespace std;

#define lambda1 1/15.0
#define lambda2 1/6.0
#define lambda3 1/9.0
#define mu1 1/6.0
#define mu2 1/5.0
#define mu3 1/4.5
#define type_client_1 0
#define type_client_2 1
#define type_client_3 2
#define e1 0
#define e2 1
#define e3 2

#define horizon 50000

typedef struct {
    int type = -1;
    double ai = 0.0;
    double ri = 0.0;
    double bi = 0.0;
    double ci = 0.0;
    double di = 0.0;
    double si = 0.0;
    double wi = 0.0;
    int numero = -1;
} client;

//variables globale
vector<client> clients_e1;
vector<client> clients_e2;
vector<client> clients_e3;

double last_ai_e1 = 0.0;
double last_ai_e2 = 0.0;
double last_ai_e3 = 0.0;

double last_ci_e1 = 0.0;
double last_ci_e2 = 0.0;
double last_ci_e3 = 0.0;

//represente le moment ou un serveur sera libre
double e1_fin = 0.0;
double e2_fin = 0.0;
double e3_fin = 0.0;


//departs
double depart_e1 = 0.0;
double depart_e2 = 0.0;
double depart_e3 = 0.0;

//cumule des durees d'attente
double cwq = 0.0;
double cwq_type1 = 0.0;
double cwq_type2 = 0.0;
double cwq_type3 = 0.0;

//cumule des durees de sejour
double cs = 0.0;
double cs_type1 = 0.0;
double cs_type2 = 0.0;
double cs_type3 = 0.0;


//ti = somme des durees ou il y'a y client de le systeme
vector<int> ti;
int pos_ti = 0.0;

MersenneTwister mt;

//fonctions
double loi_uniforme(){

    double u = 0.0;
    do{
        u = mt.random();
    } while (u < 0 || u > 1);
    return u;
}

double loi_expo(double lambda) {return ((-1.)/(lambda))*log(loi_uniforme()); }

int type_client(){
    double u = loi_uniforme();
    if (u <= 1/3.0)
        return type_client_1;
    if (u <= 2/3.0)
        return type_client_2;
    return type_client_3;
}

int choix_serveur_client_type_3(){
    double u = loi_uniforme();
    if (u <= 0.1)
        return e1;
    if (u <= 0.52)
        return e2;
    return e3;
}

bool quitter_systeme_apres_e1(){
    return loi_uniforme() <= 0.40;
}

bool quitter_systeme_apres_e2(){
    return loi_uniforme() <= 0.75;
}

void inserer(client c){
    for (int i = 0; i < clients_e3.size(); i++) {
        if(clients_e3[i].ai >= c.ai){
            clients_e3.insert(clients_e3.begin() + i, 1, c);
            return;
        }
    }
    clients_e3.push_back(c);
}

void traiter_passage_client(client c, double lambda, double mu, int ei){
    bool quitter = true;
    if(ei == e1){
        c.ri = loi_expo(lambda);
        c.ai = last_ai_e1 + c.ri;
        c.bi = fmax(last_ci_e1, c.ai);
        c.si = loi_expo(mu);
        c.di = c.bi - c.ai;
        c.ci = c.bi + c.si;
        c.wi = c.di + c.si;
        //
        cwq += c.di;
        cs += c.wi;
        if(c.type == type_client_1){
            cs_type1 += c.wi;
            cwq_type1 += c.di;
        } else{
            cs_type3 += c.wi;
            cwq_type3 += c.di;
        }
        //
        clients_e1.push_back(c);
        last_ai_e1 = c.ai;
        last_ci_e1 = c.ci;
        if(clients_e1.size() == 1){
            depart_e1 = c.ci;
        }
        else if (c.ai > depart_e1){

        }
    } else if(ei == e2){
        c.ri = loi_expo(lambda);
        c.ai = last_ai_e2 + c.ri;
        c.bi = fmax(last_ci_e2, c.ai);
        c.si = loi_expo(mu);
        c.di = c.bi - c.ai;
        c.ci = c.bi + c.si;
        c.wi = c.di + c.si;
        //
        cwq += c.di;
        cs += c.wi;
        if(c.type == type_client_2){
            cs_type2 += c.wi;
            cwq_type2 += c.di;
        } else{
            cs_type3 += c.wi;
            cwq_type3 += c.di;
        }
        //
        clients_e2.push_back(c);
        last_ai_e2 = c.ai;
        last_ci_e2 = c.ci;
    } else{
        c.ri = loi_expo(lambda);
        c.ai = last_ai_e3 + c.ri;
        c.si = loi_expo(mu);
        inserer(c);
        last_ai_e3 = c.ai;
    }

    if(c.type == type_client_1)
        quitter = quitter_systeme_apres_e1();

    if(c.type == type_client_2)
        quitter = quitter_systeme_apres_e2();

    if(!quitter){
        client c1;
        c1.si = loi_expo(mu3);
        c1.ai = c.ci;
        c1.numero = c.numero;
        c1.type = c.type;
        //cout << c1.numero << endl;
        inserer(c1);
    }
}

void run(){

    int t_client, i = 0, choix_serveur_c_t_3;
    int nbr_type1 = 0;
    int nbr_type2 = 0;
    int nbr_type3 = 0;
    /*
     * On ne considere que deux type d'evenement
     * L'arrivé et le depart
     */

    while (i < horizon)
    {
        client c;
        t_client = type_client();
        c.type = t_client;
        c.numero = i;
        if(t_client == type_client_1){
            nbr_type1++;
            traiter_passage_client(c,lambda1,mu1, e1);
        } else if (t_client == type_client_2){
            nbr_type2++;
            traiter_passage_client(c,lambda2,mu2, e2);
        } else{
            nbr_type3++;
            choix_serveur_c_t_3 = choix_serveur_client_type_3();
            if(choix_serveur_c_t_3 == e1){
                traiter_passage_client(c,lambda1,mu1, e1);
            } else if (choix_serveur_c_t_3 == e2){
                traiter_passage_client( c,lambda2,mu2, e2);
            } else{
                traiter_passage_client(c,lambda3,mu3, e3);
            }
        }
        i++;
    }

    clients_e3[0].bi = clients_e3[0].ai;
    clients_e3[0].di = 0.0;
    clients_e3[0].ci = clients_e3[0].si + clients_e3[0].bi;
    clients_e3[0].wi = clients_e3[0].si;
    for (int j = 1; j < clients_e3.size(); ++j) {
        clients_e3[j].bi = fmax(clients_e3[j].ai, clients_e3[j-1].ci);
        clients_e3[j].di = clients_e3[j].bi - clients_e3[j].ai;
        clients_e3[j].ci = clients_e3[j].bi + clients_e3[j].si;
        clients_e3[j].wi = clients_e3[j].si + clients_e3[j].di;
        //
        cwq += clients_e3[j].di;
        cs += clients_e3[j].wi;
        if(clients_e3[j].type == type_client_3){
            cwq_type3 += clients_e3[j].di;
            cs_type3 += clients_e3[j].wi;
        }
        //
        //cout << clients_e3[j].numero << "|" << clients_e3[j].ai << "|" << clients_e3[j].bi << "|" << clients_e3[j].di << "|" << clients_e3[j].ci << "|" << clients_e3[j].si << endl;
    }

    printf("---------------------------------Proportion des types de client-----------------------------------------\n");
    cout << "Type 1: " << 100*(nbr_type1/(double )horizon) << "%" << endl;
    cout << "Type 2: " << 100*(nbr_type2/(double )horizon) << "%" << endl;
    cout << "Type 3: " << 100*(nbr_type3/(double )horizon) << "%" << endl;
    printf("---------------------------------Q1-----------------------------------------\n");
    printf("Duree moyenne de sejour par type de client\n");
    printf("Type 1 : %f\n",cs_type1/nbr_type1);
    printf("Type 2 : %f\n",cs_type2/nbr_type2);
    printf("Type 3 : %f\n",cs_type3/nbr_type3);
    printf("\n---------------------------------Q2-----------------------------------------\n");
    // L = lambda * W
    printf("Lw : %f\n", (lambda1 + lambda2 + lambda3) * cs/horizon);
    printf("\n---------------------------------Q3-----------------------------------------\n");
    // L = lambda * D
    printf("Lq : %f\n", (lambda1 + lambda2 + lambda3) * cwq/horizon);
    printf("\n---------------------------------Q4-----------------------------------------\n");
    //
    cout << "Taux d'occupation (serveur 1) : "<< 100. * ((lambda1)/(mu1) + (0.1 * (lambda3)/(mu1))) << " %"  << endl;
    cout << "Taux d'occupation (serveur 2) : "<< 100. * ((lambda2)/(mu2) + (0.52 * (lambda3)/(mu2))) << " %" << endl;
    cout << "Taux d'occupation (serveur 3) : "<< 100 * ((0.38 * (lambda3)/(mu3)) + (0.6 * (lambda1)/(mu3)) + (0.25 * (lambda2)/(mu3)))<< " %"  << endl;


}


#endif //PS_RUN_H
