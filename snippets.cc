void move_warrior(int id, WarriorInfo &info) {
    Unit u = unit(id);
    Dir dir;
    if (distance_to_water(u.pos) >= u.water - 10) {
        if (info.city >= 0) {
            if (info.role == WarriorInfo::Invader or info.role == WarriorInfo::TakingOver) {
                --assigned_warriors[info.city];
            } else if (info.role == WarriorInfo::Sentinel) {
                --assigned_sentinels[info.city];
            }
            info.city = -1;
        }

        dir = nearest_water(u.pos);
    } else if (distance_to_nearest_city(u.pos) >= u.food - 10) {
        if (info.city >= 0) {
            if (info.role == WarriorInfo::Invader or info.role == WarriorInfo::TakingOver) {
                --assigned_warriors[info.city];
            } else if (info.role == WarriorInfo::Sentinel) {
                --assigned_sentinels[info.city];
            }
            info.city = -1;
        }
        dir = nearest_city(u.pos);
    } else if (info.role == WarriorInfo::Invader and info.city >= 0 and
               not city_owner_changed[info.city]) {
        if (distance_to_city(u.pos, info.city) == 0) {
            if (city_owners[info.city] == me()) {
                info.role = WarriorInfo::Sentinel;
            } else {
                info.role = WarriorInfo::TakingOver;
            }
            move_warrior(id, info);
            return;
        }
        dir = goto_city(u.pos, info.city);
    } else if (info.role == WarriorInfo::TakingOver) {
        if (city_owners[info.city] == me()) {
            --assigned_warriors[info.city];
            if (assigned_sentinels[info.city] < sentinels_per_city) {
                info.role = WarriorInfo::Sentinel;
                ++assigned_sentinels[info.city];
            } else {
                info.role = WarriorInfo::Invader;
                info.city = -1;
            }
            move_warrior(id, info);
            return;
        }
        dir = None;
        for (Dir d : random_dirs()) {
            if (city_cells[info.city].count(u.pos + d)) {
                dir = d;
                break;
            }
        }
    } else if (info.role == WarriorInfo::Sentinel) {
        if (city_owners[info.city] != me()) {
            --assigned_sentinels[info.city];
            info.role = WarriorInfo::Invader;
            info.city = -1;
            move_warrior(id, info);
            return;
        }
        dir = None;
        for (Dir d : random_dirs()) {
            if (city_cells[info.city].count(u.pos + d)) {
                dir = d;
                break;
            }
        }
    } else {
        int best_city = -1;
        int best_score = negative_infinity;
        for (unsigned int c = 0; c < city_owners.size(); ++c) {
            int score = score_city(u.pos, c);
            if (score > best_score) {
                best_score = score;
                best_city = c;
            }
        }
        info.city = best_city;
        dir = (best_city >= 0) ? goto_city(u.pos, best_city) : None;
        info.role = WarriorInfo::Invader;
    }
    if (lock(u.pos + dir)) {
        command(id, dir);
    } else if (not lock(u.pos)) {
        int d;
        for (d = 0; d < DirSize; ++d) {
            if (lock(u.pos + Dir(d))) {
                command(id, Dir(d));
                break;
            }
        }
        if (d == DirSize) {
            cerr << "Warrior " << id << " is doomed!" << endl;
        }
    }
}

Dir find_nearest_enemy(Pos starting_pos) {
    // TODO: Dijkstra
    queue<Pos> bfsq;
    bfsq.push(starting_pos);
    set<Pos> visited;
    visited.insert(starting_pos);
    Pos nearest;
    while (not bfsq.empty()) {
        Pos p = bfsq.front();
        bfsq.pop();
        Cell c = cell(p);
        if (c.type != Desert and c.type != Road) { continue; }
        if (c.id >= 0 and unit(c.id).player != me() and
            unit(c.id).type == Warrior) {
            nearest = p;
            break;
        }
        for (int d = 0; d < DirSize - 1; ++d) {
            Pos new_p = p + Dir(d);
            if (pos_ok(new_p) and visited.insert(new_p).second) {
                bfsq.push(new_p);
            }
        }
    }
    return dir_from_pos(starting_pos,
                        nearest); // A (bad) cheap approximation
}

        while (not dpq.empty()) {
            int curr_score = dpq.top().first;
            Pos curr = dpq.top().second;
            dpq.pop();
            if (visited[curr.i][curr.j]) { continue; }
            visited[curr.i][curr.j] = true;
            Pos real_pos = compute_real_pos(start_pos, origin, curr);
            if (not pos_ok(real_pos)) { continue; }
            Cell c = cell(real_pos);
            switch (c.type) {
                case Road:
                    curr_score -= 1;
                    break;
                case Desert:
                    curr_score -= 4;
                    break;
                default:
                    continue;
            }
            if (c.id != -1) {
                Unit u = unit(c.id);
                if (u.player == me() and u.type == Car) {
                    continue;
                } else if (u.player == me()) {
                    curr_score -= 100;
                } else {
                    curr_score += score_warrior(u);
                }
            }
            for (int d = 0; d < DirSize - 1; ++d) {
                Pos next = curr + Dir(d);
                if (next.i < 0 or next.j < 0 or next.i >= matrix_size or
                    next.j >= matrix_size) {
                    continue;
                }
                if (score[next.i][next.j] < curr_score) {
                    score[next.i][next.j] = curr_score;
                    prec[next.i][next.j] = curr;
                    if (curr_score >= score[current_best.i][current_best.j]) {
                        current_best = next;
                    }
                    dpq.push({score[next.i][next.j], next});
                }
            }
        }
