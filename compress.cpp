#include <bits/stdc++.h>

using namespace std;

int main() {
    string code, s;
    set<string> labels;
    while (getline(cin, s)) {
        if (s[0] == '#')
            continue;
        while (isspace(s.back()))
            s.pop_back();
        s.erase(0, find_if(s.begin(), s.end(), [](char x) {
            return !isspace(x);
        }) - s.begin());
        if (code.size())
            code.push_back(';');
        code += s;
        if (s.back() == ':' && s.size() > 7) {
            s.pop_back();
            labels.emplace(s);
        }
    }
    cerr << labels.size() << " long labels\n";
    {
        string new_code;
        for (auto x : code) {
            if (isspace(x) && new_code.size() && isspace(new_code.back()))
                continue;
            new_code.push_back(x);
        }
        code = new_code;
    }
    {
        string new_label = "AAA";
        for (auto &l : labels) {
            if (l == "main")
                continue;
            while (code.find(new_label) != string::npos) {
                for (int i = 0; i < 3; ++i) {
                    new_label[i]++;
                    if (new_label[i] > 'Z') {
                        new_label[i] = 'A';
                    } else {
                        break;
                    }
                }
            }
            string new_code;
            size_t last_i = 0;
            for (size_t i = code.find(l); i < code.size(); i = code.find(l, last_i = i + l.size())) {
                new_code += code.substr(last_i, i - last_i);
                new_code += new_label;
            }
            new_code += code.substr(last_i);
            code = new_code;
        }
    }
    const string def = "#define ";
    auto is_bad_character = [](char x) {
        return isdigit(x) || isalpha(x) || x == '_' || x == '"' || x == '$';
    };
    string preambula = "#define x(...) #__VA_ARGS__\n";
    preambula += "#define s(...) x(__VA_ARGS__)\n";
    const size_t MAX_SIZE = 99'000;
    vector<int> cs;
    for (int i = 0; i < 26 * 26 && code.size() + preambula.size() > MAX_SIZE; ++i) {
        string name = "AA";
        if (labels.count(name))
            continue;
        name[0] += i % 26;
        name[1] += i / 26;
        int best_reduction = 0;
        string best_rep;
        int c = 0;
        auto process = [&](auto &cnt) {
            for (auto &[k, v] : cnt) {
                if (find(k.begin(), k.end(), '"') != k.end())
                    continue;
                int reduction = v * (k.size() - name.size()) - def.size() - name.size() - k.size() - 2;
                if (reduction > best_reduction) {
                    best_reduction = reduction;
                    best_rep = k;
                    c = v;
                }
            }
        };
        // {
        //     unordered_map<string, int> cnt;
        //     for (auto &l : labels) {
        //         int v = 0;
        //         size_t p = 0;
        //         while (1) {
        //             p = code.find(l, p);
        //             if (p != string::npos) {
        //                 p += l.size();
        //                 ++v;
        //             } else {
        //                 break;
        //             }
        //         }
        //         --v;
        //         cnt[l] = v;
        //     }
        //     process(cnt);
        // }
        int b = 0;
        for (size_t l = 50; l > name.size() + 4; --l) {
            unordered_map<string, int> cnt;
            bool q1 = 0, q2 = 0;
            for (size_t i = 0; i < l; ++i) {
                if (code[i] == '"')
                    q2 = !q2;
                if (!q2) {
                    b += (code[i] == '(');
                    b -= (code[i] == ')');
                }
            }
            for (size_t i = 0; i + l <= code.size(); ++i) {
                if (!q1 && !q2 && b == 0 && code[i] != '(' && !isspace(code[i]) && !isspace(code[i + l - 1]) && (i == 0 || !is_bad_character(code[i - 1])) && !is_bad_character(code[i + l])) {
                    cnt[code.substr(i, l)]++;
                }
                if (code[i] == '"') {
                    q1 = !q1;
                }
                if (code[i + l] == '"') {
                    q2 = !q2;
                }
                if (!q2) {
                    b += (code[i + l] == '(');
                    b -= (code[i + l] == ')');
                }
                if (!q1) {
                    b -= (code[i] == '(');
                    b += (code[i] == ')');
                }
            }
            process(cnt);
        }
        if (!best_reduction)
            break;
        preambula += def + name + ' ' + best_rep + '\n';
        string new_code;
        size_t last_i = 0;
        for (size_t i = code.find(best_rep); i < code.size(); i = code.find(best_rep, last_i = i + best_rep.size())) {
            new_code += code.substr(last_i, i - last_i);
            if ((i == 0 || !is_bad_character(code[i - 1])) && !is_bad_character(code[i + best_rep.size()]))
                new_code += name;
            else
                new_code += best_rep;
        }
        new_code += code.substr(last_i);
        code = new_code;
        if (labels.count(best_rep)) {
            labels.erase(best_rep);
        }
        cs.push_back(c);
        cerr << best_reduction << " " << preambula.size() + code.size() << "\n";
    }
    // sort(cs.begin(), cs.end());
    // cerr << "-" << accumulate(cs.end() - 26, cs.end(), 0) << "\n";
    cerr << preambula.size() + code.size() << "\n";
    cout << preambula << "asm(s(" << code << "));\n";
}