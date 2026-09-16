#pragma once

#ifdef BHMERGER

#include <cerrno>
#include <fcntl.h>
#include <iomanip>
#include <sstream>
#include <string>
#include <unistd.h>

//! Thread- and process-safe shared output for black-hole merger events.
class BHMergerOutput {
public:
    struct ParticleState {
        long long id;
        double mass;
        double pos[3];
        double vel[3];
        double spin[3];
    };

private:
    std::string filename_;
    int precision_;
    bool enabled_;

    static bool writeLocked(const std::string& _filename,
                            const std::string& _content,
                            const bool _truncate) {
        const int flags = O_WRONLY | O_CREAT
                        | (_truncate ? O_TRUNC : O_APPEND);
        const int fd = ::open(_filename.c_str(), flags, 0666);
        if (fd < 0) return false;

        struct flock lock;
        lock.l_type = F_WRLCK;
        lock.l_whence = SEEK_SET;
        lock.l_start = 0;
        lock.l_len = 0;
        lock.l_pid = 0;

        while (::fcntl(fd, F_SETLKW, &lock) < 0) {
            if (errno != EINTR) {
                ::close(fd);
                return false;
            }
        }

        const char* data = _content.data();
        std::size_t remaining = _content.size();
        bool success = true;
        while (remaining > 0) {
            const ssize_t written = ::write(fd, data, remaining);
            if (written > 0) {
                data += written;
                remaining -= static_cast<std::size_t>(written);
            }
            else if (written < 0 && errno == EINTR) {
                continue;
            }
            else {
                success = false;
                break;
            }
        }

        lock.l_type = F_UNLCK;
        while (::fcntl(fd, F_SETLK, &lock) < 0 && errno == EINTR) {}
        ::close(fd);
        return success;
    }

    static void writeParticle(std::ostream& _out,
                              const ParticleState& _particle,
                              const int _width) {
        _out << std::setw(_width) << _particle.id
             << std::setw(_width) << _particle.mass;
        for (int k = 0; k < 3; ++k) {
            _out << std::setw(_width) << _particle.pos[k];
        }
        for (int k = 0; k < 3; ++k) {
            _out << std::setw(_width) << _particle.vel[k];
        }
        for (int k = 0; k < 3; ++k) {
            _out << std::setw(_width) << _particle.spin[k];
        }
    }

    static void writeParticleTitle(std::ostream& _out,
                                   const std::string& _prefix,
                                   const int _width) {
        _out << std::setw(_width) << (_prefix + ".id")
             << std::setw(_width) << (_prefix + ".mass")
             << std::setw(_width) << (_prefix + ".pos.x")
             << std::setw(_width) << (_prefix + ".pos.y")
             << std::setw(_width) << (_prefix + ".pos.z")
             << std::setw(_width) << (_prefix + ".vel.x")
             << std::setw(_width) << (_prefix + ".vel.y")
             << std::setw(_width) << (_prefix + ".vel.z")
             << std::setw(_width) << (_prefix + ".spin.x")
             << std::setw(_width) << (_prefix + ".spin.y")
             << std::setw(_width) << (_prefix + ".spin.z");
    }

public:
    BHMergerOutput(): filename_(), precision_(17), enabled_(false) {}

    template <class Tparticle>
    static ParticleState capture(const Tparticle& _particle) {
        ParticleState state;
        state.id = static_cast<long long>(_particle.id);
        state.mass = _particle.mass;
        state.pos[0] = _particle.pos.x;
        state.pos[1] = _particle.pos.y;
        state.pos[2] = _particle.pos.z;
        state.vel[0] = _particle.vel.x;
        state.vel[1] = _particle.vel.y;
        state.vel[2] = _particle.vel.z;
        state.spin[0] = _particle.spin.x;
        state.spin[1] = _particle.spin.y;
        state.spin[2] = _particle.spin.z;
        return state;
    }

    void open(const std::string& _filename,
              const bool _append,
              const bool _initialize_file,
              const int _width,
              const int _precision) {
        filename_ = _filename;
        precision_ = _precision;
        enabled_ = true;

        // Rank 0 alone creates a fresh file and writes its header. The caller
        // performs an MPI barrier before merger integration starts.
        if (!_append && _initialize_file) {
            std::ostringstream header;
            header << std::setw(_width) << "time";
            writeParticleTitle(header, "before1", _width);
            writeParticleTitle(header, "before2", _width);
            header << std::setw(_width) << "binary.ecc"
                   << std::setw(_width) << "binary.semi";
            writeParticleTitle(header, "remnant", _width);
            header << '\n';
            enabled_ = writeLocked(filename_, header.str(), true);
        }
    }

    bool isOpen() const {
        return enabled_;
    }

    void close() {
        enabled_ = false;
        filename_.clear();
    }

    void write(const double _time,
               const ParticleState& _before1,
               const ParticleState& _before2,
               const double _ecc,
               const double _semi,
               const ParticleState& _remnant,
               const int _width) {
        if (!enabled_) return;

#pragma omp critical(bhmerger_output)
        {
            std::ostringstream line;
            line << std::setprecision(precision_)
                 << std::setw(_width) << _time;
            writeParticle(line, _before1, _width);
            writeParticle(line, _before2, _width);
            line << std::setw(_width) << _ecc
                 << std::setw(_width) << _semi;
            writeParticle(line, _remnant, _width);
            line << '\n';
            writeLocked(filename_, line.str(), false);
        }
    }
};

#endif
