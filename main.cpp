#include "lodepng.h"
#include <bits/stdc++.h>

using namespace std;

template <class T> ostream &operator<<(ostream &o, vector<T> &arr) {
    for (auto i : arr)
        o << i << " ";
    return o;
}

template <class T> ostream &operator<<(ostream &o, vector<vector<T>> &arr) {
    for (auto &i : arr)
        o << i << " ";
    return o;
}

void compute_direction_y(const int rows, const int cols, vector<vector<int>> &g,
                         vector<vector<int>> const &mask) {
    for (int i = 0; i < cols; ++i) {
        g[0][i] = mask[0][i] == 1 ? 0 : numeric_limits<int>::max();
        for (int j = 1; j < rows; ++j) {
            if (mask[j][i] == 1) {
                g[j][i] = 0;
            } else if (g[j - 1][i] < numeric_limits<int>::max()) {
                g[j][i] = g[j - 1][i] + 1;
            }
        }
        for (int j = rows - 2; j >= 0; --j) {
            if (g[j + 1][i] < g[j][i])
                g[j][i] = g[j + 1][i] + 1;
        }
    }
}

void compute_phase_2(const int rows, const int cols,
                     vector<vector<int>> const &g, vector<vector<double>> &dt) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double best = numeric_limits<double>::max();
            for (int k = 0; k < cols; ++k) {
                int dx = abs(k - j);
                int dy = g[i][k];
                double d = dy == numeric_limits<int>::max()
                               ? numeric_limits<double>::max()
                               : (double)(dx * dx + dy * dy);
                best = min(best, d);
            }
            if (best != numeric_limits<double>::max()) {
                // cout << sqrt(best) << endl;
                dt[i][j] = sqrt(best);
            }
        }
    }
}

/// O(row*cols^2)
void compute_transform(const int rows, const int cols,
                       vector<vector<int>> const &mask,
                       vector<vector<double>> &dt) {
    vector<vector<int>> g(rows, vector<int>(cols, numeric_limits<int>::max()));
    compute_direction_y(rows, cols, g, mask);
    compute_phase_2(rows, cols, g, dt);
}

#define INF 1E20

void dt_1d(vector<double> const &f, int n, vector<double> &d) {
    vector<int> v(n);
    vector<double> z(n + 1);
    int k = 0;
    v[0] = 0;
    z[0] = -INF;
    z[1] = INF;
    for (int q = 1; q <= n - 1; q++) {
        double s =
            ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2 * q - 2 * v[k]);
        while (s <= z[k]) {
            k--;
            s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2 * q - 2 * v[k]);
        }
        k++;
        assert(k < v.size());
        v[k] = q;
        assert(k < z.size());
        z[k] = s;
        assert(k + 1 < z.size());
        z[k + 1] = INF;
    }

    k = 0;
    for (int q = 0; q <= n - 1; q++) {
        while (z[k + 1] < q)
            k++;
        d[q] = (q - v[k]) * (q - v[k]) + f[v[k]];
    }
}

// O(rows*cols)
void compute_fast(const int rows, const int cols, vector<vector<int>> const &g,
                  vector<vector<double>> &dt) {

    vector<double> f(max(rows, cols));
    assert(f.size() >= rows);
    assert(f.size() >= cols);

    // transform along columns
    for (int x = 0; x < cols; x++) {
        for (int y = 0; y < rows; y++) {
            f[y] = g[y][x] == 1 ? 0 : INF;
        }
        vector<double> d(rows);
        dt_1d(f, rows, d);
        assert(d.size() == rows);
        for (int y = 0; y < rows; y++) {
            dt[y][x] = d[y];
        }
    }

    // transform along rows
    for (int y = 0; y < rows; y++) {
        for (int x = 0; x < cols; x++) {
            double v = dt[y][x];
            f[x] = v;
        }
        vector<double> d(cols);
        dt_1d(f, cols, d);
        assert(d.size() == cols);
        for (int x = 0; x < cols; x++) {
            dt[y][x] = d[x];
        }
    }

    for (auto &i : dt)
        for (auto &j : i) {
            assert(j >= 0);
            j = sqrt(j);
        }
}

template <class T>
void save_img(string const out, const int rows, const int cols,
              vector<vector<T>> &data) {
    T vm(0);
    for (auto &i : data)
        for (auto j : i)
            if (j != numeric_limits<T>::max())
                vm = max(vm, j);

    vector<unsigned char> img_out(rows * cols * 4);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (data[i][j] == numeric_limits<T>::max()) {
                img_out[i * cols * 4 + j * 4] = 255;
                img_out[i * cols * 4 + j * 4 + 1] = 0;
                img_out[i * cols * 4 + j * 4 + 2] = 0;
                img_out[i * cols * 4 + j * 4 + 3] = 255;
            } else {
                double v = (double)data[i][j] / vm;
                v = min(v, 1.0);
                img_out[i * cols * 4 + j * 4] = 255 * v;
                img_out[i * cols * 4 + j * 4 + 1] = 255 * v;
                img_out[i * cols * 4 + j * 4 + 2] = 255 * v;
                img_out[i * cols * 4 + j * 4 + 3] = 255;
            }
        }
    }
    unsigned err = lodepng::encode(out.c_str(), img_out, cols, rows);
    if (err) {
        cout << "lodepng error: " << lodepng_error_text(err) << endl;
    }
}

void compute_gradient(const int rows, const int cols,
                      vector<vector<double>> const &dt,
                      vector<unsigned char> &out) {
    out.resize(rows * cols * 4);

    vector<vector<pair<double, double>>> dir(
        rows, vector<pair<double, double>>(cols, {0, 0}));

    double max_y(0);
    double max_x(0);
    double min_y(0);
    double min_x(0);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // derivative approximation
            if (i - 1 >= 0 && i + 1 < rows && j - 1 >= 0 && j + 1 < rows) {
                double d_y = dt[i + 1][j] - dt[i - 1][j];
                double d_x = dt[i][j + 1] - dt[i][j - 1];
                dir[i][j] = {d_y, d_x};
                max_y = max(max_y, d_y);
                max_x = max(max_x, d_x);
                min_y = min(min_y, d_y);
                min_x = min(min_x, d_x);
            }
        }
    }

    // encode direction into colour
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            auto const &[dy, dx] = dir[i][j];
            unsigned char b = 0;
            // more red for negative dy, less red for positive dy
            unsigned char r = (dy - min_y) / (max_y - min_y) * 255;
            // more green for positive dx, less green for negative dx
            unsigned char g = (dx - min_x) / (max_x - min_x) * 255;
            out[i * cols * 4 + j * 4] = r;
            out[i * cols * 4 + j * 4 + 1] = g;
            out[i * cols * 4 + j * 4 + 2] = b;
            out[i * cols * 4 + j * 4 + 3] = 255;
        }
    }
}

int main() {
    string file;
    cin >> file;
    cout << "file: " << file << endl;
    unsigned rows, cols;
    vector<unsigned char> img;

    lodepng::decode(img, cols, rows, file.c_str());

    vector<vector<int>> mask(rows, vector<int>(cols, 0));

    for (int i = 0; i < (int)rows; ++i) {
        for (int j = 0; j < (int)cols; ++j) {
            unsigned char r = img[i * cols * 4 + j * 4];
            unsigned char g = img[i * cols * 4 + j * 4 + 1];
            unsigned char b = img[i * cols * 4 + j * 4 + 2];
            if (r == 255 && g == 255 && b == 255)
                mask[i][j] = 1;
            else
                mask[i][j] = 0;
        }
    }

    save_img("imgout_mask.png", rows, cols, mask);

    vector<vector<double>> dt(
        rows, vector<double>(cols, numeric_limits<double>::max()));

    compute_transform(rows, cols, mask, dt);
    save_img("imgout.png", rows, cols, dt);

    compute_fast(rows, cols, mask, dt);
    save_img("imgout2.png", rows, cols, dt);

    vector<unsigned char> gradient;
    compute_gradient(rows, cols, dt, gradient);

    unsigned err = lodepng::encode("grad.png", gradient, cols, rows);
    if (err) {
        cout << "lodepng error: " << lodepng_error_text(err) << endl;
    }
    return 0;
}
