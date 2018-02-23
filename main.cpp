#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <stack>
#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
using namespace std;

struct str_struct {
    string a_str;
    string b_str;
};

struct coord {
    int c_i;
    int c_j;
};

struct diag_coord {
    coord coord1;
    coord coord2;
    int len;
    int sum = 0;
};

char blosum_str[24];
int blosum_matrix[24][24];

const int SIZE = 10000;
int length = 121;
//будем отбирать диагонали с длиной >= length_diag (1)
int length_diag = 2;
//отбираем наибольшие диагонали, количество которых указано в count_best_diag (2)
int count_best_diag = 2;
//сколько элементов за диагональю брать
int score = 1;
//отбираем наилучшие диагонали с учётом матрицы blosum (3)
int count_best_diag_blosum = 1;

std::ofstream ofs ("output.txt", std::ofstream::out);

void print_matrix(int** matrix, int n, int m, int n_min, int m_min) {
    for (int i = n_min; i <= n; i ++) {
        for (int j = m_min; j <= m; j ++) {
            ofs << matrix[i][j] << " ";
        }
        ofs << endl;
    }
    ofs << endl;
}

str_struct answer(const char *a, const char *b, int n, int m, int **matrix) {
    int i = n, j = m;
    str_struct str_couple;
    while (i > 0 && j > 0) {
        if (i > 0 && j > 0 && matrix[i - 1][j - 1] >= matrix[i - 1][j] && matrix[i - 1][j - 1] >= matrix[i][j - 1]) {
            i--;
            j--;
            str_couple.a_str += a[i];
            str_couple.b_str += b[j];
        } else if (i > 0 && matrix[i - 1][j] > matrix[i - 1][j - 1] && matrix[i - 1][j] > matrix[i][j - 1]) {
            i--;
            str_couple.a_str += a[i];
            str_couple.b_str += "-";
        } else {
            j--;
            str_couple.a_str += "-";
            str_couple.b_str += b[j];
        }
    }
    while (i > 0) {
        i--;
        str_couple.a_str += a[i];
        str_couple.b_str += "-";
    }
    while (j > 0) {
        j--;
        str_couple.a_str += "-";
        str_couple.b_str += b[j];
    }
    reverse(str_couple.a_str.begin(), str_couple.a_str.end());
    reverse(str_couple.b_str.begin(), str_couple.b_str.end());
    return str_couple;
}

void Nidlman_bounded(const char *a, const char *b, int **matrix) {
    int gap = 0;
    length--;
    int n = strlen(a);
    int m = strlen(b);
    for (int i = 1; i <= min(n, m+length); i++) {
        for (int j = max(i-length, 1); j <= min(i+length, m); j++) {
            int score_nid = a[i-1] == b[j-1];
            int horizontal = matrix[i][j-1] + gap;
            int vertical = matrix[i-1][j] + gap;
            int diag = matrix[i - 1][j - 1] + score_nid;
            matrix[i][j] = max(max(horizontal, vertical), diag);
        }
    }
    //print_matrix(matrix, n, m, 1, 1);
    str_struct str_couple;
    str_couple = answer(a, b, n, m, matrix);
    ofs << "Input: \n" << a  << endl << b << endl;
    ofs << endl;
    ofs << "Output: \n" << str_couple.a_str  << endl << str_couple.b_str << endl;
    ofs << endl;
}

//функция сортировки по len (по убыванию)
bool sort_arr(diag_coord i, diag_coord j) {
	return (i.len > j.len);
}

//функция сортировки по score (по убыванию)
bool sort_arr_score(diag_coord i, diag_coord j) {
	return (i.sum > j.sum);
}

int find_in_blosum (char a, char b) {
    int ind_a, ind_b, i, j;
    for (i = 0; i <= 23; i ++) {
        if (a == blosum_str[i]) {
            ind_a = i;
        }
        if (b == blosum_str[i]) {
            ind_b = i;
        }
    }
    return blosum_matrix[ind_a][ind_b];
}

void Fasta(const char *a, const char *b) {
    int n = strlen(a);
    int m = strlen(b);
    int i, j;
    int** matrix = (int**)malloc((n+1)*sizeof(int*));
    for (j = 0; j <= n; j ++) {
        matrix[j] = (int*)malloc((m+1)*sizeof(int));
    }
    for (i = 0; i <= n; i ++) {
        for (j = 0; j <= m; j ++) {
            matrix[i][j] = 0;
        }
    }
    for (i = 0; i <= n; i ++) {
        for (j = 0; j <= m; j ++) {
            if (a[i] == b[j] && (i != n) && (j != m)) {
                matrix[i][j] = 1;
            }
        }
    }
    std::vector<diag_coord> diag_array;
///------------------------------------------------------------------------------находим все диагонали, длина которых >= length_diag
    int k = 0;
    i = 0;
    //элементы над диагональю и сама диагональ
    while (1) {
        int z = 0;
        int l_diag = 0;
        coord c_start;
        coord c_end;
        diag_coord my_diag;
        int flag = 0;
        for (j = (m-k); j <= m; j ++) {
            if (matrix[z][j] == 1) {
                flag = 1;
                if ((z >= 1) && (j >= 1) && matrix[z-1][j-1] == 0) {
                    l_diag = 0;
                }
                if (l_diag == 0) {
                    c_start.c_i = z;
                    c_start.c_j = j;
                }
                l_diag ++;
                if (l_diag >= length_diag) {
                    my_diag.coord1 = c_start;
                }
            }
            z++;
            if (flag && l_diag >= length_diag && ((z > n && j > m) || (matrix[z][j+1] == 0))) {
                flag = 0;
                c_end.c_i = z-1;
                c_end.c_j = j;
                my_diag.coord2 = c_end;
                my_diag.len = l_diag;
                diag_array.push_back(my_diag);
            }
        }
        i = 0;
        k ++;
        if (z == j) {
            //прошли диагональ
            break;
        }
    }
    //элементы под диагональю
    for (i = 1; i <= n; i++) {
        int z = i;
        j = 0;
        int l_diag = 0;
        coord c_start;
        coord c_end;
        diag_coord my_diag;
        int flag = 0;
        while ((z <= n) && (j <= m)) {
            if (matrix[z][j] == 1) {
                flag = 1;
                if ((z >= 1) && (j >= 1) && matrix[z-1][j-1] == 0) {
                    l_diag = 0;
                }
                if (l_diag == 0) {
                    c_start.c_i = z;
                    c_start.c_j = j;
                }
                l_diag ++;
                if (l_diag >= length_diag) {
                    my_diag.coord1 = c_start;
                }
            }
            z ++;
            j ++;
            if (flag && l_diag >= length_diag && ((z > n && j > m) || (matrix[z][j] == 0))) {
                flag = 0;
                c_end.c_i = z-1;
                c_end.c_j = j-1;
                my_diag.coord2 = c_end;
                my_diag.len = l_diag;
                diag_array.push_back(my_diag);
            }
        }
    }
///------------------------------------------------------------------------------обходим вектор с диагоналями и ищем лучшие диагонали
    vector<diag_coord>::iterator it_diag;
    sort(diag_array.begin(), diag_array.end(), sort_arr);   ///отсортировали вектор с диагоналями по убыванию
    ///оставляем только лучшие диагонали в векторе; их количество = count_best_diag
    while (diag_array.size() != count_best_diag) {
        diag_array.pop_back();
    }
    std::vector<diag_coord> score_diag_array;
    ///оставили только те диагонали, которые можем взять с учётом score
    for (it_diag = diag_array.begin(); it_diag != diag_array.end(); it_diag++) {
        //cout << "len:" << (*it_diag).len << ", i_start: " << (*it_diag).coord1.c_i << ", j_start: " << (*it_diag).coord1.c_j << "; i_end: "  << (*it_diag).coord2.c_i << ", j_end: " << (*it_diag).coord2.c_j << endl;
        if (((*it_diag).coord1.c_i - score) >= 0 && ((*it_diag).coord1.c_j - score) >= 0 && ((*it_diag).coord2.c_i + score) <= n && ((*it_diag).coord2.c_j + score) <= m) {
            (*it_diag).coord1.c_i -= score;
            (*it_diag).coord1.c_j -= score;
            (*it_diag).coord2.c_i += score;
            (*it_diag).coord2.c_j += score;
            score_diag_array.push_back(*it_diag);
        }
    }
    ///подсчитываем сумму для каждой из диагоналей
    for (it_diag = score_diag_array.begin(); it_diag != score_diag_array.end(); it_diag++) {
        int j = (*it_diag).coord1.c_j;
        for (i = (*it_diag).coord1.c_i; i <= (*it_diag).coord2.c_i; i++) {
            /*cout << a[i] << " ";
            cout << b[j] << " ";
            cout << find_in_blosum(a[i], b[i]) << endl;*/
            (*it_diag).sum += find_in_blosum(a[i], b[i]);
            j ++;
        }
    }

    sort(score_diag_array.begin(), score_diag_array.end(), sort_arr_score);   ///отсортировали вектор с диагоналями по сумме score
    ///оставляем только лучшие диагонали в векторе; их количество = count_best_diag_blosum
    while (score_diag_array.size() != count_best_diag_blosum) {
        score_diag_array.pop_back();
    }
    int i_max = 0;
    int j_max = 0;
    for (it_diag = score_diag_array.begin(); it_diag != score_diag_array.end(); it_diag++) {
        if ((*it_diag).coord1.c_i > i_max) {
            i_max = (*it_diag).coord1.c_i;
        }
        if ((*it_diag).coord1.c_j > i_max) {
            j_max = (*it_diag).coord1.c_j;
        }
    }
    int** n_matrix = (int**)malloc((n+1)*sizeof(int*));
    for (int j = 0; j <= n; j ++) {
        n_matrix[j] = (int*)malloc((m+1)*sizeof(int));
    }
    for (int i = 0; i <= n; i ++) {
        for (int j = 0; j <= m; j ++) {
            n_matrix[i][j] = 0;
        }
    }
    length = max(i_max, j_max);
    Nidlman_bounded(a, b, n_matrix);
    return;
}

void input_FASTA(const char a[SIZE], const char b[SIZE]) {
    int n, m;
    if (strlen(a) > strlen(b)) {
        n = strlen(a);
        m = strlen(b);
    } else {
        m = strlen(a);
        n = strlen(b);
    }
    if (strlen(a) > strlen(b)) {
        Fasta(a, b);
    } else {
        Fasta(b, a);
    }
}

int main() {
    FILE *in=fopen("blosum.txt", "rt");
    for (int i = 0; i <= 23; i++) {
        fscanf(in, "%c", &blosum_str[i]);
    }
    for (int i = 0; i <= 23; i++) {
        for (int j = 0; j <= 23; j++) {
            fscanf(in, "%d", &blosum_matrix[i][j]);
        }
    }
    int count_n;
    cout << "Enter the number of sequences" << endl;
    scanf("%d", &count_n);
    freopen("input.txt", "a+", stdin);
    string arr_of_str[count_n];
    for (int i = 0; i < count_n; i ++) {
        char a[SIZE];
        gets(a);
        arr_of_str[i] = a;
    }
    for (int i = 0; i < count_n; i ++) {
        for (int j = i+1; j < count_n; j ++) {
            input_FASTA(arr_of_str[i].c_str(), arr_of_str[j].c_str());
        }
    }
    return 0;
}
