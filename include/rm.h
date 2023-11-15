#ifndef EXP_RM_H
#define EXP_RM_H

#include <random>

#include "IMs.h"
#include "OPIM_new.h"
#include "Heap.h"
#include "rm_rrset.h"


double calc_bound_RM(Graph &G, int64 T, RMRRContainer &RRI, std::vector<double> &q_R) {
    double sum = 0;
    for (int i = 0; i < G.n; i++) {
        double value = 0;
        for (int j = 0; j < T; ++j) {
            double value_v = 0;
            for (auto rr: RRI.covered[j * G.n + i]) {
                value_v += q_R[rr];
            }
            value = std::max(value, value_v);
        }
        sum += value;
    }
    return sum;
}

void rounding_RM(Graph &G, std::vector<std::vector<bi_node>> &bases, std::vector<int64> &frac_x,
                 std::vector<double> &q_R, RMRRContainer &RRI, int64 t_max, std::vector<bi_node> &seeds) {
    std::vector<int64> pos(G.n);
    std::vector<int64> pos1(G.n);
    for (auto &i: bases[0]) {
        //first: round second:node
        pos[i.second] = i.first;
    }
    for (int i = 0; i < bases.size() - 1; ++i) {
        for (auto &j: bases[i + 1]) {
            //first: round second:node
            pos1[j.second] = j.first;
        }
        for (int j = 0; j < pos.size(); ++j) {
            //j: node
            int x = pos[j], y = pos1[j];
            //x, y: round
            if (x != y) {
                // round (x,j) and (y,j)
                double dx = 0, dy = 0;
                for (auto rr: RRI.covered[G.n * x + j]) {
                    dx += q_R[rr];
                }
                dx *= (double) t_max / (double) (t_max - frac_x[G.n * x + j]);
                for (auto rr: RRI.covered[G.n * y + j]) {
                    dy += q_R[rr];
                }
                dy *= (double) t_max / (double) (t_max - frac_x[G.n * y + j]);
                int x_ = x, y_ = y, swap_flag = 1;
                if (dx < dy) { //|| (dx == dy && random_real() > (double) (i + 1) / (i + 2))
                    std::swap(x_, y_);
                    pos[j] = y;
                    swap_flag = i + 1;
                } else {
                    pos1[j] = x;
                }
                for (auto rr: RRI.covered[G.n * y_ + j]) {
                    q_R[rr] *= (double) (t_max - frac_x[G.n * y_ + j] + swap_flag) /
                               (double) (t_max - frac_x[G.n * y_ + j]);
                }
                frac_x[G.n * y_ + j] -= swap_flag;
                for (auto rr: RRI.covered[G.n * x_ + j]) {
                    q_R[rr] *= (double) (t_max - frac_x[G.n * x_ + j] - swap_flag) /
                               (double) (t_max - frac_x[G.n * x_ + j]);
                }
                frac_x[G.n * x_ + j] += swap_flag;
                //end of rounding x and y
            }
        }
    }
    for (int i = 0; i < pos.size(); ++i) {
        // first : node second : round
        seeds.emplace_back(i, pos[i]);
    }
}

double CGreedy_RM(Graph &G, RMRRContainer &RRI, int64 T, int64 t_max, std::vector<bi_node> &bi_seeds) {

    //the fractional solution (use integers to avoid float error)
    std::vector<int64> frac_x(G.n * T, 0);
    //temporary varible
    double Fx = 0, tight_bound = RRI.numOfRRsets();
    std::vector<double> q_R(RRI.numOfRRsets(), 1);

    //priority queue for lazy sampling
    HeapArray<CMax<double, std::pair<std::pair<unsigned, unsigned>, unsigned>>> Q{};
    Q.val = new double[G.n * T];
    Q.ids = new std::pair<std::pair<unsigned, unsigned>, unsigned>[G.n * T];
    Q.k = 0;

    std::vector<std::vector<bi_node>> bases(t_max);

    double cur = clock();
    for (int t = 1; t <= t_max; t++) {
        tight_bound = std::min(tight_bound, Fx + calc_bound_RM(G, T, RRI, q_R));
        Q.k = 0;
        double sum = 0;
        for (int j = 0; j < G.n; j++) {
            double value_max = 0;
            for (int i = 0; i < T; i++) {
                double value_v = 0;
                for (auto rr: RRI.covered[G.n * i + j]) {
                    value_v += q_R[rr];
                }
                value_v *= (double) t_max / (double) (t_max - frac_x[G.n * i + j]);
                value_max = std::max(value_max, value_v);
                Q.k++;
                heap_push<CMax<double, std::pair<std::pair<unsigned, unsigned>, unsigned>>>(Q.k, Q.val, Q.ids,
                                                                                            value_v, std::make_pair(
                                std::make_pair(i, j), 0));
            }
            sum += value_max;
        }
        tight_bound = std::min(tight_bound, Fx + sum);
        std::vector<bool> node_selected(G.n, false);
        for (int i = 0; i < G.n; i++) {
            while (Q.k != 0) {
                //j: round l:node
                int64 j = Q.ids[0].first.first;
                int64 l = Q.ids[0].first.second;
                int64 it_round = Q.ids[0].second;
                heap_pop<CMax<double, std::pair<std::pair<unsigned, unsigned>, unsigned>>>(Q.k, Q.val, Q.ids);
                Q.k -= 1;
                if (node_selected[l]) continue;
                if (it_round == i) {
                    //choose
                    node_selected[l] = true;
                    for (auto rr: RRI.covered[G.n * j + l]) {
                        double q_R_old = q_R[rr];
                        q_R[rr] *= (double) (t_max - frac_x[G.n * j + l] - 1) /
                                   (double) (t_max - frac_x[G.n * j + l]);
                        Fx += q_R_old - q_R[rr];
                    }
                    frac_x[G.n * j + l] += 1;
                    bases[t - 1].emplace_back(j, l);
                    break;
                } else {
                    double value_v = 0;
                    for (auto rr: RRI.covered[G.n * j + l]) {
                        value_v += q_R[rr];
                    }
                    value_v *= (double) t_max / (double) (t_max - frac_x[G.n * j + l]);
                    Q.k++;
                    heap_push<CMax<double, std::pair<std::pair<unsigned, unsigned>, unsigned>>>(Q.k, Q.val, Q.ids,
                                                                                                value_v, std::make_pair(
                                    std::make_pair(j, l), i));
                }
            }
        }
    }
    std::cout << "alg time = " << (clock() - cur) / CLOCKS_PER_SEC;

    tight_bound = std::min(tight_bound, Fx + calc_bound_RM(G, T, RRI, q_R));

    rounding_RM(G, bases, frac_x, q_R, RRI, t_max, bi_seeds);
    delete[] Q.val;
    delete[] Q.ids;
    return tight_bound;
}

double CGreedy_RM_PM(Graph &G, RMRRContainer &RRI, int64 T, int64 t_max, std::vector<bi_node> &bi_seeds) {

    //the fractional solution (use integers to avoid float error)
    std::vector<int64> frac_x(G.n * T, 0);
    //temporary varible
    double Fx = 0;
    std::vector<double> q_R(RRI.numOfRRsets(), 1);

    std::vector<std::vector<bi_node>> bases(t_max);

    double cur = clock();
    for (int t = 1; t <= t_max; t++) {
        for (int i = 0; i < G.n; i++) {
            double value = -1;
            int selected_round = -1;
            for (int j = 0; j < T; j++) {
                double value_v = 0;
                for (auto rr: RRI.covered[G.n * j + i]) {
                    value_v += q_R[rr];
                }
                value_v *= (double) t_max / (double) (t_max - frac_x[G.n * j + i]);
                if (value_v > value) selected_round = j, value = value_v;
            }
            for (auto rr: RRI.covered[G.n * selected_round + i]) {
                double q_R_old = q_R[rr];
                q_R[rr] *= (double) (t_max - frac_x[G.n * selected_round + i] - 1) /
                           (double) (t_max - frac_x[G.n * selected_round + i]);
                Fx += q_R_old - q_R[rr];
            }
            frac_x[G.n * selected_round + i] += 1;
            bases[t - 1].emplace_back(selected_round, i);
        }
    }
    std::cout << "alg time = " << (clock() - cur) / CLOCKS_PER_SEC;

    double tight_bound = Fx + calc_bound_RM(G, T, RRI, q_R);

    rounding_RM(G, bases, frac_x, q_R, RRI, t_max, bi_seeds);
    return tight_bound;
}


#endif //EXP_RM_H
