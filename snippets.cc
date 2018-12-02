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
