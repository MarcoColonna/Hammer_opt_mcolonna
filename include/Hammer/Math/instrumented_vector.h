#pragma once

#include <cstdint>
#include <initializer_list>
#include <type_traits>
#include <vector>
#include <limits>
#include <algorithm>
#include <utility>
#include <cstdio>
#include <unistd.h>

/** @brief instrumented version of std::vector that keeps track of the maximum
 * size and capacity.
 *
 * It's enough to instrument the destructor and the methods that can
 * potentially reduce the size or capacity (so we don't miss the maximum).
 */
template <class T, class ALLOC = std::allocator<T>, class BASE = std::vector<T, ALLOC>>
class instrumented_vector : public BASE {
public:
    // duplicate all the standard typedefs
    using size_type = typename BASE::size_type;
    using difference_type = typename BASE::difference_type;
    using value_type = typename BASE::value_type;
    using pointer = typename BASE::pointer;
    using const_pointer = typename BASE::const_pointer;
    using reference = typename BASE::reference;
    using const_reference = typename BASE::const_reference;
    using iterator = typename BASE::iterator;
    using const_iterator = typename BASE::const_iterator;
    using reverse_iterator = typename BASE::reverse_iterator;
    using const_reverse_iterator = typename BASE::const_reverse_iterator;
    using allocator_type = typename BASE::allocator_type;

private:
    // keep track of max size and capacity of instances during their lifetime
    size_type m_max_size = std::numeric_limits<size_type>::min();
    size_type m_max_cap = std::numeric_limits<size_type>::min();
    // keep track of min/max of the above for all instances of the class
    inline static size_type s_min_max_size = std::numeric_limits<size_type>::max();
    inline static size_type s_max_max_size = std::numeric_limits<size_type>::min();
    inline static size_type s_min_max_cap = std::numeric_limits<size_type>::max();
    inline static size_type s_max_max_cap = std::numeric_limits<size_type>::min();
    constexpr static std::size_t histsz = 32;
    inline static std::size_t s_hist_size[histsz] = {};
    inline static std::size_t s_hist_cap[histsz] = {};
    
    /// keep track of per-instance max size and capacity
    void update_members() noexcept
    {
        if (BASE::size() > m_max_size) m_max_size = BASE::size();
        if (BASE::capacity() > m_max_cap) m_max_cap = BASE::capacity();
    }
    /// keep track of min/max of per-class max size and capacity
    void update_static_members() noexcept
    {
        if (m_max_size > s_max_max_size) s_max_max_size = m_max_size;
        if (m_max_size < s_min_max_size) s_min_max_size = m_max_size;
        if (m_max_cap > s_max_max_cap) s_max_max_cap = m_max_cap;
        if (m_max_cap < s_min_max_cap) s_min_max_cap = m_max_cap;
        ++s_hist_size[std::min(m_max_size, histsz - 1)];
        ++s_hist_cap[std::min(m_max_cap, histsz - 1)];
    }
    static void print_at_exit(void)
    {
        const int pid = getpid();
        constexpr auto bufsz = 256;
        char fname[bufsz];
        int n = std::snprintf(fname, bufsz, "instrumented_vector.%u.log", pid);
        fname[std::min(n, bufsz - 1)] = 0;
        FILE *file = std::fopen(fname, "a");
        std::fprintf(file,
                     "[DEBUG] In %s, pid %u: size min %zu max %zu capacity "
                     "min %zu "
                     "max %zu\n",
                     __PRETTY_FUNCTION__, pid, std::size_t(s_min_max_size),
                     std::size_t(s_max_max_size), std::size_t(s_min_max_cap),
                     std::size_t(s_max_max_cap));
        constexpr char legend[4] = { ' ', '/', '\\', 'X' };
        std::fprintf(file,
                     "[DEBUG] In %s: histograms (legend: '%c' size, '%c' "
                     "capacity, '%c' both)\n",
                     __PRETTY_FUNCTION__, legend[1], legend[2], legend[3]);
        const auto max_sz = std::max(std::size_t(1), *std::max_element(std::begin(s_hist_size), std::end(s_hist_size)));
        const auto max_cap = std::max(std::size_t(1), *std::max_element(std::begin(s_hist_cap), std::end(s_hist_cap)));
        constexpr auto termsz = 80;
        for (std::size_t i = 0; histsz != i; ++i) {
            unsigned col0 = std::fprintf(file, "sz %9zu cap %9zu|%s%2zu|",
                                         s_hist_size[i], s_hist_cap[i],
                                         (i == histsz - 1 ? ">=" : "  "), i);
            const std::size_t linelen = termsz - col0;
            for (unsigned col = 0; linelen != col; ++col) {
                fputc(legend[1 * (col ==
                                  std::min((128 * linelen * s_hist_size[i] +
                                            64) / (128 * max_sz),
                                           linelen - 1)) +
                             2 * (col ==
                                  std::min((128 * linelen * s_hist_cap[i] +
                                            64) / (128 * max_cap),
                                           linelen - 1))],
                      file);
            }
            std::fputc('\n', file);
        }
        std::fclose(file);
    }

    class at_exit {
    private:
        void (*m_fn)(void);
    public:
        at_exit(void (*fn)(void)) : m_fn(fn) {}
        ~at_exit() { m_fn(); }
    };

    static inline at_exit s_atexit = {print_at_exit};

public:
    // some constructors need special treatment...
    instrumented_vector(const instrumented_vector&) = default;
    instrumented_vector(instrumented_vector&& other)
            : BASE(std::move(other)), m_max_size(other.m_max_size),
              m_max_cap(other.m_max_cap)
    {
        other.m_max_size = other.m_max_cap = 0;
    }
    instrumented_vector(const BASE& other) : BASE(other) {}
    instrumented_vector(BASE&& other) : BASE(std::move(other)) {}
    template <class U, class = std::enable_if_t<std::is_same<
                               std::initializer_list<value_type>, U>::value>>
    instrumented_vector(U&& ilist,
                        const allocator_type& alloc = allocator_type())
            : BASE(std::forward<U>(ilist), alloc)
    {}
    // can reuse the other constructors verbatim
    using BASE::BASE;

    // instrument destructor and operations that can shrink size
    ~instrumented_vector() {
        update_members();
        update_static_members();
        // we have to appear to be using s_atexit in some code that's
        // actually used, or the compiler will optimize it away (the
        // destructor does work fine for that)
        static_cast<void>(&s_atexit);
    }

    instrumented_vector& operator=(const instrumented_vector& other) {
        update_members();
        BASE::operator=(other);
        return *this;
    }

    void swap(instrumented_vector& other) noexcept(noexcept(BASE::swap(std::declval<BASE&>())))
    {
        update_members();
        std::swap(m_max_size, other.m_max_size);
        std::swap(m_max_cap, other.m_max_cap);
        BASE::swap(other);
        update_members();
    }

    instrumented_vector& operator=(instrumented_vector&& other) {
        if (&other != this) instrumented_vector::swap(other);
        return *this;
    }

    void clear()
    {
        update_members();
        BASE::clear();
    }

    template <class... ARGS>
    void assign(ARGS&&... args)
    {
        update_members();
        BASE::assign(std::forward<ARGS>(args)...);
    }

    template <class... ARGS>
    void resize(ARGS&&... args)
    {
        update_members();
        BASE::resize(std::forward<ARGS>(args)...);
    }

    template <class... ARGS>
    iterator erase(ARGS&&... args)
    {
        update_members();
        return BASE::erase(std::forward<ARGS>(args)...);
    }
    
    void shrink_to_fit()
    {
        update_members();
        BASE::shrink_to_fit();
    }
};
