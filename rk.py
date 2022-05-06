from numpy import zeros


# Metoda Eulera - układ trzech równań
def Euler_method_3(X, Y, Z, X_0, Y_0, Z_0, N_t, h):
    czas = [0]
    for i in range(int((N_t) * (1 / h)) - 1):
        czas.append(czas[i] + h)
    a, b, c = zeros(len(czas)), zeros(len(czas)), zeros(len(czas))
    a[0], b[0], c[0] = X_0, Y_0, Z_0
    for j in range(len(czas) - 1):
        a[j + 1] = a[j] + X(czas[j], a[j], b[j], c[j]) * h
        b[j + 1] = b[j] + Y(czas[j], a[j], b[j], c[j]) * h
        c[j + 1] = c[j] + Z(czas[j], a[j], b[j], c[j]) * h
    return czas, a, b, c


# Metoda Rungego-Kutty 4. rzędu - układ trzech równań
def rk4_method_3(X, Y, Z, X_0, Y_0, Z_0, N_t, h):
    czas = [0]
    for i in range(int((N_t) * (1 / h)) - 1):
        czas.append(czas[i] + h)
    a, b, c = zeros(len(czas)), zeros(len(czas)), zeros(len(czas))
    a[0], b[0], c[0] = X_0, Y_0, Z_0
    for j in range(len(czas) - 1):
        k_1 = X(czas[j], a[j], b[j], c[j]) * h
        l_1 = Y(czas[j], a[j], b[j], c[j]) * h
        m_1 = Z(czas[j], a[j], b[j], c[j]) * h

        k_2 = X(czas[j] + 0.5 * h, a[j] + 0.5 * k_1, b[j] + 0.5 * l_1, c[j] + 0.5 * m_1) * h
        l_2 = Y(czas[j] + 0.5 * h, a[j] + 0.5 * k_1, b[j] + 0.5 * l_1, c[j] + 0.5 * m_1) * h
        m_2 = Z(czas[j] + 0.5 * h, a[j] + 0.5 * k_1, b[j] + 0.5 * l_1, c[j] + 0.5 * m_1) * h

        k_3 = X(czas[j] + 0.5 * h, a[j] + 0.5 * k_2, b[j] + 0.5 * l_2, c[j] + 0.5 * m_2) * h
        l_3 = Y(czas[j] + 0.5 * h, a[j] + 0.5 * k_2, b[j] + 0.5 * l_2, c[j] + 0.5 * m_2) * h
        m_3 = Z(czas[j] + 0.5 * h, a[j] + 0.5 * k_2, b[j] + 0.5 * l_2, c[j] + 0.5 * m_2) * h

        k_4 = X(czas[j] + h, a[j] + k_3, b[j] + l_3, c[j] + m_3) * h
        l_4 = Y(czas[j] + h, a[j] + k_3, b[j] + l_3, c[j] + m_3) * h
        m_4 = Z(czas[j] + h, a[j] + k_3, b[j] + l_3, c[j] + m_3) * h

        a[j + 1] = a[j] + (1 / 6) * (k_1 + 2 * k_2 + 2 * k_3 + k_4)
        b[j + 1] = b[j] + (1 / 6) * (l_1 + 2 * l_2 + 2 * l_3 + l_4)
        c[j + 1] = c[j] + (1 / 6) * (m_1 + 2 * m_2 + 2 * m_3 + m_4)
    return czas, a, b, c


# Metoda Rungego-Kutty 4. rzędu - układ czterech równań
def rk4_method_4(X, Y, Z, U, X_0, Y_0, Z_0, U_0, N_t, h):
    czas = [0]
    for i in range(int((N_t) * (1 / h)) - 1):
        czas.append(czas[i] + h)
    a, b, c, d = zeros(len(czas)), zeros(len(czas)), zeros(len(czas)), zeros(len(czas))
    a[0], b[0], c[0], d[0] = X_0, Y_0, Z_0, U_0
    for j in range(len(czas) - 1):
        k_1 = X(czas[j], a[j], b[j], c[j], d[j]) * h
        l_1 = Y(czas[j], a[j], b[j], c[j], d[j]) * h
        m_1 = Z(czas[j], a[j], b[j], c[j], d[j]) * h
        n_1 = U(czas[j], a[j], b[j], c[j], d[j]) * h

        k_2 = X(czas[j] + 0.5 * h, a[j] + 0.5 * k_1, b[j] + 0.5 * l_1, c[j] + 0.5 * m_1, d[j] + 0.5 * n_1) * h
        l_2 = Y(czas[j] + 0.5 * h, a[j] + 0.5 * k_1, b[j] + 0.5 * l_1, c[j] + 0.5 * m_1, d[j] + 0.5 * n_1) * h
        m_2 = Z(czas[j] + 0.5 * h, a[j] + 0.5 * k_1, b[j] + 0.5 * l_1, c[j] + 0.5 * m_1, d[j] + 0.5 * n_1) * h
        n_2 = U(czas[j] + 0.5 * h, a[j] + 0.5 * k_1, b[j] + 0.5 * l_1, c[j] + 0.5 * m_1, d[j] + 0.5 * n_1) * h

        k_3 = X(czas[j] + 0.5 * h, a[j] + 0.5 * k_2, b[j] + 0.5 * l_2, c[j] + 0.5 * m_2, d[j] + 0.5 * n_2) * h
        l_3 = Y(czas[j] + 0.5 * h, a[j] + 0.5 * k_2, b[j] + 0.5 * l_2, c[j] + 0.5 * m_2, d[j] + 0.5 * n_2) * h
        m_3 = Z(czas[j] + 0.5 * h, a[j] + 0.5 * k_2, b[j] + 0.5 * l_2, c[j] + 0.5 * m_2, d[j] + 0.5 * n_2) * h
        n_3 = U(czas[j] + 0.5 * h, a[j] + 0.5 * k_2, b[j] + 0.5 * l_2, c[j] + 0.5 * m_2, d[j] + 0.5 * n_2) * h

        k_4 = X(czas[j] + h, a[j] + k_3, b[j] + l_3, c[j] + m_3, d[j] + n_3) * h
        l_4 = Y(czas[j] + h, a[j] + k_3, b[j] + l_3, c[j] + m_3, d[j] + n_3) * h
        m_4 = Z(czas[j] + h, a[j] + k_3, b[j] + l_3, c[j] + m_3, d[j] + n_3) * h
        n_4 = U(czas[j] + h, a[j] + k_3, b[j] + l_3, c[j] + m_3, d[j] + n_3) * h

        a[j + 1] = a[j] + (1 / 6) * (k_1 + 2 * k_2 + 2 * k_3 + k_4)
        b[j + 1] = b[j] + (1 / 6) * (l_1 + 2 * l_2 + 2 * l_3 + l_4)
        c[j + 1] = c[j] + (1 / 6) * (m_1 + 2 * m_2 + 2 * m_3 + m_4)
        d[j + 1] = d[j] + (1 / 6) * (n_1 + 2 * n_2 + 2 * n_3 + n_4)
    return czas, a, b, c, d
