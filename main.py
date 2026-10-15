import sys
import math
import numpy as np

# Configuración de codificación UTF-8 para consola en Windows
if sys.platform == 'win32':
    try:
        if sys.stdout is not None:
            sys.stdout.reconfigure(encoding='utf-8')
        if sys.stderr is not None:
            sys.stderr.reconfigure(encoding='utf-8')
    except Exception:
        pass

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QGridLayout, QLabel, QPushButton, QComboBox, QDoubleSpinBox,
    QGroupBox
)
from PyQt5.QtCore import Qt, QTimer, QPointF, QRectF
from PyQt5.QtGui import QPainter, QPen, QColor, QFont, QFontMetrics, QBrush
import pyqtgraph as pg


# ==============================================================================
# 1. MOTOR FÍSICO Y ANALÍTICO DE TRANSITORIOS (v2.0 Riguroso)
# ==============================================================================
class TransientEngine:
    """
    Motor analítico exacto para circuitos RLC Serie, RLC Paralelo, RC y RL.
    Soporta 4 estados de fuente: Directa (+Vs), Invertida (-Vs), Corto (0V) y Abierta (Homogénea).
    Conserva rigurosamente la continuidad física: v_C(0+) = v_C(0-) e i_L(0+) = i_L(0-).
    """
    SW_DIRECT = 0    # +Vs (Forzada directa)
    SW_INVERT = 1    # -Vs (Forzada invertida)
    SW_SHORT = 2     # 0V  (Homogénea con fuente en corto)
    SW_OPEN = 3      # Abierta (Homogénea pura sin rama de fuente)

    def __init__(self):
        self.topology = 'RLC_SERIE'  # 'RLC_SERIE', 'RLC_PARALELO', 'RC_SERIE', 'RL_SERIE'

        # Parámetros base
        self.Vs = 12.0       # Magnitud de fuente principal (V)
        self.R1 = 100.0      # Resistencia principal / acoplamiento (Ohm)
        self.R2 = 220.0      # Resistencia paralela / descarga (Ohm)
        self.L = 0.05        # Inductancia (H) -> 50 mH
        self.C = 10e-6       # Capacitancia (F) -> 10 uF

        # Estado de llaves
        self.sw1_state = self.SW_DIRECT
        self.sw2_closed = True

        # Variables de estado físicas continuas
        self.vC = 0.0
        self.iL = 0.0
        self.vL = 0.0
        self.v_node = 0.0

        # Corrientes de rama exactas según leyes de Kirchhoff
        self.i_source = 0.0
        self.iR1 = 0.0
        self.iR2 = 0.0
        self.i_load = 0.0
        self.iC = 0.0

        # Equivalente activo visto por la red reactiva
        self.v_th = 12.0
        self.r_eq = 100.0
        self.is_homogeneous = False

        # Tiempos
        self.t_sim = 0.0
        self.t_local = 0.0

        # Condiciones iniciales capturadas en t0+
        self.vC_0 = 0.0
        self.iL_0 = 0.0
        self.dvC_0 = 0.0
        self.diL_0 = 0.0

        # Parámetros característicos
        self.vC_inf = 0.0
        self.iL_inf = 0.0
        self.tau = 0.001
        self.alpha = 0.0
        self.w0 = 0.0
        self.wd = 0.0
        self.zeta = 0.0
        self.regime_type = ""
        self.roots_str = ""

        # Coeficientes de evaluación y tokens de ecuaciones
        self._vC_coeffs = None
        self._iL_coeffs = None
        self.eq_symbolic_tokens = []
        self.eq_numeric_tokens = []
        self.eq_tau_tokens = []
        self.eq_sol_v_tokens = []
        self.eq_sol_i_tokens = []

        self.recalculate_transient()

    def get_source_voltage_value(self):
        if self.sw1_state == self.SW_DIRECT:
            return self.Vs
        elif self.sw1_state == self.SW_INVERT:
            return -self.Vs
        else:
            return 0.0

    def recalculate_transient(self):
        """
        Captura continuidad vC(0+) = vC(0-) e iL(0+) = iL(0-) y resuelve
        analíticamente la EDO (Forzada u Homogénea según S1).
        """
        self.t_local = 0.0
        self.vC_0 = self.vC
        self.iL_0 = self.iL

        L = max(self.L, 1e-9)
        C = max(self.C, 1e-12)
        R1 = max(self.R1, 1e-3)
        R2 = max(self.R2, 1e-3)

        self.is_homogeneous = (self.sw1_state in (self.SW_SHORT, self.SW_OPEN))
        Vs_val = self.get_source_voltage_value()

        # ------------------------------------------------------------------
        # DETERMINACIÓN DEL CIRCUITO EQUIVALENTE SEGÚN TOPOLOGÍA Y LLAVES
        # ------------------------------------------------------------------
        if self.topology in ('RLC_SERIE', 'RC_SERIE', 'RL_SERIE'):
            if self.sw1_state == self.SW_OPEN:
                self.v_th = 0.0
                if not self.sw2_closed:
                    self.r_eq = float('inf')
                else:
                    self.r_eq = R2
            elif self.sw1_state == self.SW_SHORT:
                self.v_th = 0.0
                if self.sw2_closed:
                    self.r_eq = (R1 * R2) / (R1 + R2)
                else:
                    self.r_eq = R1
            else:
                if not self.sw2_closed:
                    self.v_th = Vs_val
                    self.r_eq = R1
                else:
                    self.v_th = Vs_val * (R2 / (R1 + R2))
                    self.r_eq = (R1 * R2) / (R1 + R2)

        elif self.topology == 'RLC_PARALELO':
            if self.sw1_state == self.SW_OPEN:
                self.v_th = 0.0
                self.r_eq = R2 if self.sw2_closed else float('inf')
            elif self.sw1_state == self.SW_SHORT:
                self.v_th = 0.0
                self.r_eq = (R1 * R2) / (R1 + R2) if self.sw2_closed else R1
            else:
                self.v_th = Vs_val
                self.r_eq = (R1 * R2) / (R1 + R2) if self.sw2_closed else R1

        R = self.r_eq
        Vth = self.v_th

        # ------------------------------------------------------------------
        # CASO 1: RC SERIE (1° ORDEN)
        # ------------------------------------------------------------------
        if self.topology == 'RC_SERIE':
            self.iL = 0.0
            self.iL_0 = 0.0
            if math.isinf(R):
                self.tau = 1.0
                self.alpha = 0.0
                self.w0 = 0.0
                self.zeta = 1.0
                self.vC_inf = self.vC_0
                self.dvC_0 = 0.0
                self.regime_type = "RC Serie — Circuito Abierto (Carga Retenida)"
                self.roots_str = "s = 0"
                self._vC_coeffs = ('CONST', self.vC_0)
                self._iL_coeffs = ('CONST', 0.0)
                self._build_tokens_rc_open()
                self.step(0.0)
                return

            self.tau = R * C
            self.alpha = 1.0 / self.tau
            self.w0 = 0.0
            self.zeta = 1.0
            self.vC_inf = Vth
            self.iL_inf = 0.0
            self.dvC_0 = (Vth - self.vC_0) / (R * C)
            self.regime_type = "RC Serie — Homogénea (Descarga)" if self.is_homogeneous else "RC Serie — Forzada (Carga)"
            self.roots_str = f"s = -{self.alpha:.1f} s⁻¹"

            A = self.vC_0 - Vth
            I0 = (Vth - self.vC_0) / R
            self._vC_coeffs = ('1ST', Vth, A, self.tau)
            self._iL_coeffs = ('1ST', 0.0, I0, self.tau)
            self._build_tokens_1st_order('RC', R, L, C, Vth, A, I0)

        # ------------------------------------------------------------------
        # CASO 2: RL SERIE (1° ORDEN)
        # ------------------------------------------------------------------
        elif self.topology == 'RL_SERIE':
            self.vC = 0.0
            self.vC_0 = 0.0
            if math.isinf(R):
                self.tau = 1.0
                self.alpha = float('inf')
                self.w0 = 0.0
                self.zeta = 1.0
                self.iL_inf = 0.0
                self.vC_inf = 0.0
                self.diL_0 = 0.0
                self.regime_type = "RL Serie — Circuito Abierto (Desconectado)"
                self.roots_str = "i_L = 0"
                self._vC_coeffs = ('CONST', 0.0)
                self._iL_coeffs = ('CONST', 0.0)
                self.step(0.0)
                return

            self.tau = L / R
            self.alpha = 1.0 / self.tau
            self.w0 = 0.0
            self.zeta = 1.0
            self.iL_inf = Vth / R
            self.vC_inf = 0.0
            self.diL_0 = (Vth - R * self.iL_0) / L
            self.regime_type = "RL Serie — Homogénea (Descarga)" if self.is_homogeneous else "RL Serie — Forzada (Carga)"
            self.roots_str = f"s = -{self.alpha:.1f} s⁻¹"

            A_i = self.iL_0 - self.iL_inf
            V_L0 = Vth - R * self.iL_0
            self._iL_coeffs = ('1ST', self.iL_inf, A_i, self.tau)
            self._vC_coeffs = ('1ST', 0.0, V_L0, self.tau)
            self._build_tokens_1st_order('RL', R, L, C, Vth, V_L0, A_i)

        # ------------------------------------------------------------------
        # CASO 3: RLC SERIE (2° ORDEN)
        # ------------------------------------------------------------------
        elif self.topology == 'RLC_SERIE':
            if math.isinf(R):
                self.alpha = 0.0
                self.w0 = 1.0 / math.sqrt(L * C)
                self.tau = 1.0
                self.zeta = float('inf')
                self.vC_inf = self.vC_0
                self.iL_inf = 0.0
                self.regime_type = "RLC Serie — Circuito Abierto (Lazo Cortado)"
                self.roots_str = "i_L = 0"
                self._vC_coeffs = ('CONST', self.vC_0)
                self._iL_coeffs = ('CONST', 0.0)
                self.step(0.0)
                return

            self.alpha = R / (2.0 * L)
            self.w0 = 1.0 / math.sqrt(L * C)
            self.tau = 1.0 / self.alpha if self.alpha > 1e-9 else (2.0 * math.pi / self.w0)
            self.zeta = self.alpha / self.w0

            self.vC_inf = Vth
            self.iL_inf = 0.0
            self.dvC_0 = self.iL_0 / C
            self.diL_0 = (Vth - R * self.iL_0 - self.vC_0) / L

            self._vC_coeffs, self.eq_sol_v_tokens = self._solve_2nd_order(
                self.vC_0, self.dvC_0, self.vC_inf, self.alpha, self.w0, "v", "C", "V"
            )
            self._iL_coeffs, self.eq_sol_i_tokens = self._solve_2nd_order(
                self.iL_0, self.diL_0, self.iL_inf, self.alpha, self.w0, "i", "L", "A"
            )
            self._build_tokens_rlc_serie(R, L, C, Vth)

        # ------------------------------------------------------------------
        # CASO 4: RLC PARALELO (2° ORDEN)
        # ------------------------------------------------------------------
        elif self.topology == 'RLC_PARALELO':
            self.w0 = 1.0 / math.sqrt(L * C)
            if math.isinf(R):
                self.alpha = 0.0
                self.tau = (2.0 * math.pi) / self.w0
            else:
                self.alpha = 1.0 / (2.0 * R * C)
                self.tau = 1.0 / self.alpha
            self.zeta = self.alpha / self.w0

            In = 0.0 if self.is_homogeneous else (Vs_val / R1)
            self.vC_inf = 0.0
            self.iL_inf = In

            i_R0 = 0.0 if math.isinf(R) else (self.vC_0 / R)
            self.dvC_0 = (In - i_R0 - self.iL_0) / C
            self.diL_0 = self.vC_0 / L

            self._vC_coeffs, self.eq_sol_v_tokens = self._solve_2nd_order(
                self.vC_0, self.dvC_0, self.vC_inf, self.alpha, self.w0, "v", "C", "V"
            )
            self._iL_coeffs, self.eq_sol_i_tokens = self._solve_2nd_order(
                self.iL_0, self.diL_0, self.iL_inf, self.alpha, self.w0, "i", "L", "A"
            )
            self._build_tokens_rlc_paralelo(R, L, C, In)

        self.step(0.0)

    def _solve_2nd_order(self, x0, dx0, x_inf, alpha, w0, var_base, var_sub, unit):
        e0 = x0 - x_inf
        ratio = alpha / w0 if w0 > 0 else 1.0
        prec = ".2f" if unit == "V" else ".3f"
        c_var = "#00E676" if unit == "V" else "#FFEA00"

        tokens = [('sub', var_base, var_sub, c_var), ('t', "(t) = ", "#E2E8F0")]
        if abs(x_inf) > 1e-5:
            tokens.append(('t', f"{x_inf:{prec}} + ", "#38BDF8"))

        if alpha < 1e-6:
            self.regime_type = "RLC — Oscilador LC Puro (α = 0, ζ = 0)"
            self.wd = w0
            A1 = e0
            A2 = dx0 / w0
            coeffs = ('UNDER', x_inf, A1, A2, 0.0, w0)
            self.roots_str = f"s₁,₂ = ± j{w0:.1f} rad/s"
            tokens.extend([
                ('t', f"[ ({A1:+{prec}})·cos({w0:.1f}t) + ({A2:+{prec}})·sin({w0:.1f}t) ]", c_var),
                ('t', f"  [{unit}]", "#94A3B8")
            ])
        elif abs(ratio - 1.0) < 1e-3:
            self.regime_type = "RLC — Críticamente Amortiguado (ζ = 1)"
            A1 = e0
            A2 = dx0 + alpha * e0
            coeffs = ('CRIT', x_inf, A1, A2, alpha)
            self.roots_str = f"s₁ = s₂ = -{alpha:.1f} s⁻¹"
            tokens.extend([
                ('t', f"[ ({A1:+{prec}}) + ({A2:+.2e})·t ] · ", c_var),
                ('exp', f"-{alpha:.1f}t", "#FF9100"),
                ('t', f"  [{unit}]", "#94A3B8")
            ])
        elif alpha > w0:
            self.regime_type = "RLC — Sobreamortiguado (ζ > 1)"
            disc = math.sqrt(alpha**2 - w0**2)
            s1 = -alpha + disc
            s2 = -alpha - disc
            A1 = (dx0 - s2 * e0) / (s1 - s2)
            A2 = e0 - A1
            coeffs = ('OVER', x_inf, A1, A2, s1, s2)
            self.roots_str = f"s₁ = {s1:.1f} s⁻¹, s₂ = {s2:.1f} s⁻¹"
            tokens.extend([
                ('t', f"({A1:+{prec}})·", c_var),
                ('exp', f"{s1:.1f}t", "#FF9100"),
                ('t', f" + ({A2:+{prec}})·", c_var),
                ('exp', f"{s2:.1f}t", "#FF9100"),
                ('t', f"  [{unit}]", "#94A3B8")
            ])
        else:
            self.regime_type = "RLC — Subamortiguado / Oscilatorio (ζ < 1)"
            wd = math.sqrt(w0**2 - alpha**2)
            self.wd = wd
            A1 = e0
            A2 = (dx0 + alpha * e0) / wd
            coeffs = ('UNDER', x_inf, A1, A2, alpha, wd)
            self.roots_str = f"s₁,₂ = -{alpha:.1f} ± j{wd:.1f} rad/s"
            tokens.extend([
                ('exp', f"-{alpha:.1f}t", "#FF9100"),
                ('t', f" · [ ({A1:+{prec}})·cos({wd:.1f}t) + ({A2:+{prec}})·sin({wd:.1f}t) ]", c_var),
                ('t', f"  [{unit}]", "#94A3B8")
            ])

        return coeffs, tokens

    def _build_tokens_rc_open(self):
        self.eq_symbolic_tokens = [
            ('frac', [('t', 'dv', '#00E676'), ('sub', '', 'C', '#00E676')], [('t', 'dt', '#E2E8F0')], '#E2E8F0'),
            ('t', ' = 0   [Circuito Abierto: i(t) = 0]', '#00E5FF')
        ]
        self.eq_numeric_tokens = [
            ('frac', [('t', 'dv', '#00E676'), ('sub', '', 'C', '#00E676')], [('t', 'dt', '#E2E8F0')], '#E2E8F0'),
            ('t', ' = 0.00 V/s', '#00E5FF')
        ]
        self.eq_tau_tokens = [
            ('t', 'τ = ∞   |   Carga retenida constante en capacitor (Sin camino de descarga)   |   s = 0', '#FFEA00')
        ]
        self.eq_sol_v_tokens = [
            ('sub', 'v', 'C', '#00E676'),
            ('t', f'(t) = {self.vC_0:+.2f} V  (Carga Retenida Constante)', '#00E676')
        ]
        self.eq_sol_i_tokens = [
            ('t', 'i(t) = 0.000 A', '#FFEA00')
        ]

    def _build_tokens_1st_order(self, kind, R, L, C, Vth, A_v, A_i):
        rhs_color = "#00E5FF" if self.is_homogeneous else "#FF9100"
        if kind == 'RC':
            self.eq_symbolic_tokens = [
                ('frac', [('t', 'dv', '#00E676'), ('sub', '', 'C', '#00E676')], [('t', 'dt', '#E2E8F0')], '#E2E8F0'),
                ('t', ' + ', '#E2E8F0'),
                ('frac', [('t', '1', '#E2E8F0')], [('sub', 'R', 'eq', '#FF9100'), ('t', '·C', '#00E676')], '#E2E8F0'),
                ('t', ' · ', '#E2E8F0'),
                ('sub', 'v', 'C', '#00E676'),
                ('t', '(t) = ', '#E2E8F0'),
                ('t', '0' if self.is_homogeneous else 'V_th / (R_eq·C)', rhs_color)
            ]
            rhs_num = "0" if self.is_homogeneous else f"{Vth/(R*C):.2e}"
            self.eq_numeric_tokens = [
                ('frac', [('t', 'dv', '#00E676'), ('sub', '', 'C', '#00E676')], [('t', 'dt', '#E2E8F0')], '#E2E8F0'),
                ('t', f' + {1.0/self.tau:.1f} · ', '#E2E8F0'),
                ('sub', 'v', 'C', '#00E676'),
                ('t', '(t) = ', '#E2E8F0'),
                ('t', rhs_num, rhs_color)
            ]
            self.eq_tau_tokens = [
                ('t', f'τ = {self.format_eng(self.tau, "s")}  (5τ = {self.format_eng(5*self.tau, "s")})   │   {self.roots_str}', '#E2E8F0')
            ]
        else:
            self.eq_symbolic_tokens = [
                ('frac', [('t', 'di', '#FFEA00'), ('sub', '', 'L', '#FFEA00')], [('t', 'dt', '#E2E8F0')], '#E2E8F0'),
                ('t', ' + ', '#E2E8F0'),
                ('frac', [('sub', 'R', 'eq', '#FF9100')], [('t', 'L', '#FFEA00')], '#E2E8F0'),
                ('t', ' · ', '#E2E8F0'),
                ('sub', 'i', 'L', '#FFEA00'),
                ('t', '(t) = ', '#E2E8F0'),
                ('t', '0' if self.is_homogeneous else 'V_th / L', rhs_color)
            ]
            rhs_num = "0" if self.is_homogeneous else f"{Vth/L:.2e}"
            self.eq_numeric_tokens = [
                ('frac', [('t', 'di', '#FFEA00'), ('sub', '', 'L', '#FFEA00')], [('t', 'dt', '#E2E8F0')], '#E2E8F0'),
                ('t', f' + {1.0/self.tau:.1f} · ', '#E2E8F0'),
                ('sub', 'i', 'L', '#FFEA00'),
                ('t', '(t) = ', '#E2E8F0'),
                ('t', rhs_num, rhs_color)
            ]
            self.eq_tau_tokens = [
                ('t', f'τ = {self.format_eng(self.tau, "s")}  (5τ = {self.format_eng(5*self.tau, "s")})   │   {self.roots_str}', '#E2E8F0')
            ]

        self.eq_sol_v_tokens = [
            ('sub', 'v', 'C' if kind == 'RC' else 'L', '#00E676'),
            ('t', f'(t) = {self.vC_inf:.2f} + ({A_v:+.2f}) · ', '#00E676'),
            ('exp', f'-t / {self.format_eng(self.tau, "s")}', '#FF9100'),
            ('t', '  [V]', '#94A3B8')
        ]
        self.eq_sol_i_tokens = [
            ('sub', 'i', 'L' if kind == 'RL' else 'C', '#FFEA00'),
            ('t', f'(t) = {self.iL_inf:.3f} + ({A_i:+.3f}) · ', '#FFEA00'),
            ('exp', f'-t / {self.format_eng(self.tau, "s")}', '#FF9100'),
            ('t', '  [A]', '#94A3B8')
        ]

    def _build_tokens_rlc_serie(self, R, L, C, Vth):
        rhs_color = "#00E5FF" if self.is_homogeneous else "#FF9100"
        self.eq_symbolic_tokens = [
            ('frac', [('sup', 'd', '2', '#E2E8F0'), ('sub', 'v', 'C', '#00E676')], [('t', 'dt', '#E2E8F0'), ('sup', '', '2', '#E2E8F0')], '#E2E8F0'),
            ('t', ' + ', '#E2E8F0'),
            ('frac', [('sub', 'R', 'eq', '#FF9100')], [('t', 'L', '#FFEA00')], '#E2E8F0'),
            ('t', ' · ', '#E2E8F0'),
            ('frac', [('t', 'dv', '#00E676'), ('sub', '', 'C', '#00E676')], [('t', 'dt', '#E2E8F0')], '#E2E8F0'),
            ('t', ' + ', '#E2E8F0'),
            ('frac', [('t', '1', '#E2E8F0')], [('t', 'L·C', '#00E5FF')], '#E2E8F0'),
            ('t', ' · ', '#E2E8F0'),
            ('sub', 'v', 'C', '#00E676'),
            ('t', '(t) = ', '#E2E8F0'),
            ('t', '0' if self.is_homogeneous else 'V_th / (L·C)', rhs_color)
        ]
        rhs_val = "0" if self.is_homogeneous else f"{(self.w0**2)*Vth:.2e}"
        self.eq_numeric_tokens = [
            ('frac', [('sup', 'd', '2', '#E2E8F0'), ('sub', 'v', 'C', '#00E676')], [('t', 'dt', '#E2E8F0'), ('sup', '', '2', '#E2E8F0')], '#E2E8F0'),
            ('t', f' + {2*self.alpha:.1f} · ', '#E2E8F0'),
            ('frac', [('t', 'dv', '#00E676'), ('sub', '', 'C', '#00E676')], [('t', 'dt', '#E2E8F0')], '#E2E8F0'),
            ('t', f' + {self.w0**2:.2e} · ', '#E2E8F0'),
            ('sub', 'v', 'C', '#00E676'),
            ('t', '(t) = ', '#E2E8F0'),
            ('t', rhs_val, rhs_color)
        ]
        self.eq_tau_tokens = [
            ('t', f'τ = {self.format_eng(self.tau, "s")}  (5τ = {self.format_eng(5*self.tau, "s")})   │   '
                  f'α = {self.alpha:.1f} s⁻¹   │   ω₀ = {self.w0:.1f} rad/s (ζ = {self.zeta:.2f})   │   {self.roots_str}', '#E2E8F0')
        ]

    def _build_tokens_rlc_paralelo(self, R, L, C, In):
        rhs_color = "#00E5FF" if self.is_homogeneous else "#FF9100"
        self.eq_symbolic_tokens = [
            ('frac', [('sup', 'd', '2', '#E2E8F0'), ('sub', 'i', 'L', '#FFEA00')], [('t', 'dt', '#E2E8F0'), ('sup', '', '2', '#E2E8F0')], '#E2E8F0'),
            ('t', ' + ', '#E2E8F0'),
            ('frac', [('t', '1', '#E2E8F0')], [('sub', 'R', 'p', '#FF9100'), ('t', '·C', '#00E676')], '#E2E8F0'),
            ('t', ' · ', '#E2E8F0'),
            ('frac', [('t', 'di', '#FFEA00'), ('sub', '', 'L', '#FFEA00')], [('t', 'dt', '#E2E8F0')], '#E2E8F0'),
            ('t', ' + ', '#E2E8F0'),
            ('frac', [('t', '1', '#E2E8F0')], [('t', 'L·C', '#00E5FF')], '#E2E8F0'),
            ('t', ' · ', '#E2E8F0'),
            ('sub', 'i', 'L', '#FFEA00'),
            ('t', '(t) = ', '#E2E8F0'),
            ('t', '0' if self.is_homogeneous else 'I_N / (L·C)', rhs_color)
        ]
        rhs_val = "0" if self.is_homogeneous else f"{(self.w0**2)*In:.2e}"
        self.eq_numeric_tokens = [
            ('frac', [('sup', 'd', '2', '#E2E8F0'), ('sub', 'i', 'L', '#FFEA00')], [('t', 'dt', '#E2E8F0'), ('sup', '', '2', '#E2E8F0')], '#E2E8F0'),
            ('t', f' + {2*self.alpha:.1f} · ', '#E2E8F0'),
            ('frac', [('t', 'di', '#FFEA00'), ('sub', '', 'L', '#FFEA00')], [('t', 'dt', '#E2E8F0')], '#E2E8F0'),
            ('t', f' + {self.w0**2:.2e} · ', '#E2E8F0'),
            ('sub', 'i', 'L', '#FFEA00'),
            ('t', '(t) = ', '#E2E8F0'),
            ('t', rhs_val, rhs_color)
        ]
        self.eq_tau_tokens = [
            ('t', f'τ = {self.format_eng(self.tau, "s")}  (5τ = {self.format_eng(5*self.tau, "s")})   │   '
                  f'α = {self.alpha:.1f} s⁻¹   │   ω₀ = {self.w0:.1f} rad/s (ζ = {self.zeta:.2f})   │   {self.roots_str}', '#E2E8F0')
        ]

    def _eval_analytic(self, coeffs, t):
        if coeffs is None:
            return 0.0
        kind = coeffs[0]
        if kind == 'CONST':
            return coeffs[1]
        elif kind == '1ST':
            _, x_inf, A, tau = coeffs
            decay = math.exp(-min(t / max(tau, 1e-12), 100.0))
            return x_inf + A * decay
        elif kind == 'CRIT':
            _, x_inf, A1, A2, alpha = coeffs
            decay = math.exp(-min(alpha * t, 100.0))
            return x_inf + (A1 + A2 * t) * decay
        elif kind == 'OVER':
            _, x_inf, A1, A2, s1, s2 = coeffs
            d1 = math.exp(max(s1 * t, -100.0))
            d2 = math.exp(max(s2 * t, -100.0))
            return x_inf + A1 * d1 + A2 * d2
        elif kind == 'UNDER':
            _, x_inf, A1, A2, alpha, wd = coeffs
            decay = math.exp(-min(alpha * t, 100.0))
            return x_inf + decay * (A1 * math.cos(wd * t) + A2 * math.sin(wd * t))
        return 0.0

    def step(self, dt):
        """
        Avanza el reloj y calcula todas las corrientes de rama con conservación
        estricta de las leyes de Kirchhoff.
        """
        self.t_sim += dt
        self.t_local += dt

        val_v = self._eval_analytic(self._vC_coeffs, self.t_local)
        val_i = self._eval_analytic(self._iL_coeffs, self.t_local)

        Vs_val = self.get_source_voltage_value()

        if self.topology == 'RLC_SERIE':
            self.vC = val_v
            self.iL = val_i
            self.i_load = val_i
            self.iC = val_i

            # Tensión en el nodo N1 (divisor con R2 y carga)
            if self.sw1_state == self.SW_OPEN:
                self.i_source = 0.0
                if self.sw2_closed:
                    self.v_node = -self.R2 * self.i_load
                    self.iR2 = -self.i_load
                else:
                    self.v_node = 0.0
                    self.iR2 = 0.0
            else:
                self.v_node = self.v_th - self.r_eq * self.i_load
                self.iR2 = (self.v_node / self.R2) if self.sw2_closed else 0.0
                self.i_source = self.iR2 + self.i_load

            self.vL = self.v_node - self.vC

        elif self.topology == 'RC_SERIE':
            self.vC = val_v
            self.iL = 0.0
            self.v_node = self.vC

            if self.sw2_closed:
                self.iR2 = self.v_node / self.R2
            else:
                self.iR2 = 0.0

            if self.sw1_state == self.SW_OPEN:
                self.i_source = 0.0
            elif self.sw1_state == self.SW_SHORT:
                self.i_source = -self.v_node / self.R1
            else:
                self.i_source = (Vs_val - self.v_node) / self.R1

            # Por KCL en nodo N1: i_C = i_source - i_R2
            self.iC = self.i_source - self.iR2
            self.i_load = self.iC
            self.vL = 0.0

        elif self.topology == 'RL_SERIE':
            self.vC = 0.0
            self.iL = val_i
            self.i_load = val_i
            self.iC = 0.0

            if self.sw1_state == self.SW_OPEN:
                self.i_source = 0.0
                if self.sw2_closed:
                    self.v_node = -self.R2 * self.iL
                    self.iR2 = -self.iL
                else:
                    self.v_node = 0.0
                    self.iR2 = 0.0
            else:
                self.v_node = self.v_th - self.r_eq * self.iL
                self.iR2 = (self.v_node / self.R2) if self.sw2_closed else 0.0
                self.i_source = self.iR2 + self.iL

            self.vL = self.v_node

        elif self.topology == 'RLC_PARALELO':
            self.vC = val_v
            self.iL = val_i
            self.v_node = self.vC
            self.vL = self.vC

            if self.sw1_state == self.SW_OPEN:
                self.i_source = 0.0
            elif self.sw1_state == self.SW_SHORT:
                self.i_source = -self.vC / self.R1
            else:
                self.i_source = (Vs_val - self.vC) / self.R1

            self.iR2 = (self.vC / self.R2) if self.sw2_closed else 0.0
            # KCL: i_source = i_R2 + i_L + i_C
            self.iC = self.i_source - self.iR2 - self.iL
            self.i_load = self.iL + self.iC

        self.iR1 = self.i_source

    def force_steady_state(self):
        self.t_local = 25.0 * max(self.tau, 1e-6)
        self.vC = self.vC_inf
        self.iL = self.iL_inf
        self.step(0.0)

    def is_in_steady_state(self):
        if self.alpha < 1e-6:
            return False
        return self.t_local >= 5.0 * max(self.tau, 1e-9)

    @staticmethod
    def format_eng(val, unit=""):
        if math.isinf(val):
            return f"∞ {unit}"
        if val == 0:
            return f"0 {unit}"
        abs_v = abs(val)
        if abs_v >= 1.0:
            return f"{val:.3f} {unit}"
        elif abs_v >= 1e-3:
            return f"{val*1e3:.2f} m{unit}"
        elif abs_v >= 1e-6:
            return f"{val*1e6:.2f} µ{unit}"
        else:
            return f"{val*1e9:.2f} n{unit}"


# ==============================================================================
# 2. WIDGET DE ECUACIONES MATEMÁTICAS EN GRAN FORMATO (FRACCIONES REALES)
# ==============================================================================
class MathEquationWidget(QWidget):
    """
    Renderiza la Ecuación Diferencial Simbólica, Numérica, Constantes y Soluciones
    con fracciones apiladas verticales sin colisiones, tipografía matemática Cambria Math,
    y badge luminoso sin superposiciones.
    """
    def __init__(self, engine: TransientEngine):
        super().__init__()
        self.engine = engine
        self.setMinimumHeight(255)

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        painter.setRenderHint(QPainter.TextAntialiasing)

        w = self.width()
        h = self.height()

        # Fondo pizarra matemática
        painter.fillRect(0, 0, w, h, QColor("#0D111A"))
        painter.setPen(QPen(QColor("#243044"), 1.5))
        painter.drawRoundedRect(1, 1, w - 2, h - 2, 8, 8)

        # -------------------------------------------------------------
        # CABECERA: TÍTULO, CONDICIONES INICIALES Y BADGE (SIN COLISIÓN)
        # -------------------------------------------------------------
        badge_w = 245
        badge_h = 24
        badge_rect = QRectF(w - badge_w - 14, 7, badge_w, badge_h)

        if self.engine.is_homogeneous:
            painter.setBrush(QBrush(QColor("#064E3B")))
            painter.setPen(QPen(QColor("#10B981"), 1.5))
            badge_txt = "EDO HOMOGÉNEA ( = 0 )"
            badge_col = QColor("#6EE7B7")
        else:
            painter.setBrush(QBrush(QColor("#451A03")))
            painter.setPen(QPen(QColor("#F59E0B"), 1.5))
            badge_txt = "EDO FORZADA (NO HOMOGÉNEA)"
            badge_col = QColor("#FDE68A")

        painter.drawRoundedRect(badge_rect, 5, 5)
        painter.setFont(QFont("Segoe UI", 8, QFont.Bold))
        painter.setPen(badge_col)
        painter.drawText(badge_rect, Qt.AlignCenter, badge_txt)

        # Título a la izquierda
        painter.setFont(QFont("Segoe UI", 10, QFont.Bold))
        painter.setPen(QColor("#38BDF8"))
        fm_title = QFontMetrics(painter.font())
        title_txt = self.engine.regime_type
        title_w = fm_title.horizontalAdvance(title_txt)
        painter.drawText(16, 24, title_txt)

        # Condiciones iniciales compactas con verificación de colisión estricta
        painter.setFont(QFont("Segoe UI", 8, QFont.Bold))
        painter.setPen(QColor("#94A3B8"))
        fm_ic = QFontMetrics(painter.font())
        ic_txt = f"CI: [ v_C(0⁺) = {self.engine.vC_0:+.2f} V ,  i_L(0⁺) = {self.engine.iL_0*1e3:+.1f} mA ]"
        ic_w = fm_ic.horizontalAdvance(ic_txt)
        badge_left = w - badge_w - 14

        # Colocar CI únicamente si queda espacio entre el título y el badge con margen seguro
        if 16 + title_w + 20 + ic_w < badge_left:
            painter.drawText(int(badge_left - ic_w - 16), 24, ic_txt)
        elif 16 + title_w + 12 + ic_w < badge_left:
            painter.drawText(int(16 + title_w + 12), 24, ic_txt)

        # Línea divisoria
        painter.setPen(QPen(QColor("#1E293B"), 1))
        painter.drawLine(14, 38, w - 14, 38)

        # Columna de inicio común para todos los tokens
        col_x = 118

        # -------------------------------------------------------------
        # FILA 1: EDO SIMBÓLICA (Con fracciones verticales)
        # -------------------------------------------------------------
        y_sym = 66
        painter.setFont(QFont("Segoe UI", 8, QFont.Bold))
        painter.setPen(QColor("#64748B"))
        painter.drawText(16, y_sym + 4, "SIMBÓLICA:")
        self._draw_tokens(painter, col_x, y_sym, self.engine.eq_symbolic_tokens, base_pt=13)

        # -------------------------------------------------------------
        # FILA 2: EDO NUMÉRICA REAL (Con fracciones verticales)
        # -------------------------------------------------------------
        y_num = 112
        painter.setFont(QFont("Segoe UI", 8, QFont.Bold))
        painter.setPen(QColor("#64748B"))
        painter.drawText(16, y_num + 4, "NUMÉRICA:")
        self._draw_tokens(painter, col_x, y_num, self.engine.eq_numeric_tokens, base_pt=13)

        # Línea divisoria
        painter.setPen(QPen(QColor("#1E293B"), 1))
        painter.drawLine(14, 138, w - 14, 138)

        # -------------------------------------------------------------
        # FILA 3: CONSTANTE DE TIEMPO τ Y RAÍCES
        # -------------------------------------------------------------
        y_tau = 160
        painter.setFont(QFont("Segoe UI", 8, QFont.Bold))
        painter.setPen(QColor("#64748B"))
        painter.drawText(16, y_tau + 4, "PARÁMETROS:")
        self._draw_tokens(painter, col_x, y_tau, self.engine.eq_tau_tokens, base_pt=9.5)

        # -------------------------------------------------------------
        # FILA 4: SOLUCIÓN TENSIÓN v(t)
        # -------------------------------------------------------------
        y_v = 192
        painter.setFont(QFont("Segoe UI", 8, QFont.Bold))
        painter.setPen(QColor("#00E676"))
        painter.drawText(16, y_v + 4, "TENSIÓN:")
        self._draw_tokens(painter, col_x, y_v, self.engine.eq_sol_v_tokens, base_pt=11.5)

        # -------------------------------------------------------------
        # FILA 5: SOLUCIÓN CORRIENTE i(t)
        # -------------------------------------------------------------
        y_i = 224
        painter.setFont(QFont("Segoe UI", 8, QFont.Bold))
        painter.setPen(QColor("#FFEA00"))
        painter.drawText(16, y_i + 4, "CORRIENTE:")
        self._draw_tokens(painter, col_x, y_i, self.engine.eq_sol_i_tokens, base_pt=11.5)

    def _get_math_font(self, pt, italic=False, bold=True):
        f = QFont("Cambria Math", int(pt))
        f.setBold(bold)
        f.setItalic(italic)
        return f

    def _measure_subtokens(self, painter, toks, pt):
        w = 0
        for tok in toks:
            kind = tok[0]
            if kind == 't':
                painter.setFont(self._get_math_font(pt))
                w += QFontMetrics(painter.font()).horizontalAdvance(tok[1])
            elif kind in ('sub', 'sup'):
                painter.setFont(self._get_math_font(pt, italic=True))
                w_b = QFontMetrics(painter.font()).horizontalAdvance(tok[1]) if tok[1] else 0
                painter.setFont(self._get_math_font(max(pt - 3, 8)))
                w_s = QFontMetrics(painter.font()).horizontalAdvance(tok[2])
                w += w_b + w_s + 2
        return w

    def _draw_subtokens(self, painter, x, baseline, toks, pt, is_num=True):
        cur_x = x
        for tok in toks:
            kind = tok[0]
            if kind == 't':
                painter.setFont(self._get_math_font(pt))
                painter.setPen(QColor(tok[2]))
                painter.drawText(int(cur_x), int(baseline), tok[1])
                cur_x += QFontMetrics(painter.font()).horizontalAdvance(tok[1])
            elif kind == 'sub':
                painter.setPen(QColor(tok[3]))
                if tok[1]:
                    painter.setFont(self._get_math_font(pt, italic=True))
                    painter.drawText(int(cur_x), int(baseline), tok[1])
                    cur_x += QFontMetrics(painter.font()).horizontalAdvance(tok[1])
                painter.setFont(self._get_math_font(max(pt - 3, 8)))
                sub_y = baseline + 1 if is_num else baseline + 3
                painter.drawText(int(cur_x), int(sub_y), tok[2])
                cur_x += QFontMetrics(painter.font()).horizontalAdvance(tok[2]) + 2
            elif kind == 'sup':
                painter.setPen(QColor(tok[3]))
                if tok[1]:
                    painter.setFont(self._get_math_font(pt, italic=True))
                    painter.drawText(int(cur_x), int(baseline), tok[1])
                    cur_x += QFontMetrics(painter.font()).horizontalAdvance(tok[1])
                painter.setFont(self._get_math_font(max(pt - 3, 8)))
                sup_y = baseline - 5 if is_num else baseline - 4
                painter.drawText(int(cur_x), int(sup_y), tok[2])
                cur_x += QFontMetrics(painter.font()).horizontalAdvance(tok[2]) + 2
        return cur_x

    def _draw_frac(self, painter, cur_x, bar_y, num_toks, den_toks, col, pt):
        frac_pt = max(pt - 2, 10)
        w_num = self._measure_subtokens(painter, num_toks, frac_pt)
        w_den = self._measure_subtokens(painter, den_toks, frac_pt)
        bar_w = max(w_num, w_den) + 8
        painter.setPen(QPen(QColor(col), 1.5))
        painter.drawLine(int(cur_x + 1), int(bar_y), int(cur_x + 1 + bar_w), int(bar_y))
        num_x = cur_x + 1 + (bar_w - w_num) / 2.0
        self._draw_subtokens(painter, num_x, bar_y - 4, num_toks, frac_pt, is_num=True)
        den_x = cur_x + 1 + (bar_w - w_den) / 2.0
        self._draw_subtokens(painter, den_x, bar_y + 12, den_toks, frac_pt, is_num=False)
        return cur_x + bar_w + 6

    def _draw_tokens(self, painter, x, y_center, tokens, base_pt=13):
        cur_x = x
        for tok in tokens:
            kind = tok[0]
            if kind == 't':
                _, text, col = tok
                painter.setFont(self._get_math_font(base_pt))
                painter.setPen(QColor(col))
                painter.drawText(int(cur_x), int(y_center + 4), text)
                cur_x += QFontMetrics(painter.font()).horizontalAdvance(text)

            elif kind == 'sub':
                _, base, sub, col = tok
                painter.setPen(QColor(col))
                if base:
                    painter.setFont(self._get_math_font(base_pt, italic=True))
                    painter.drawText(int(cur_x), int(y_center + 4), base)
                    cur_x += QFontMetrics(painter.font()).horizontalAdvance(base)
                painter.setFont(self._get_math_font(max(base_pt - 3, 8)))
                painter.drawText(int(cur_x), int(y_center + 8), sub)
                cur_x += QFontMetrics(painter.font()).horizontalAdvance(sub) + 2

            elif kind == 'sup':
                _, base, sup, col = tok
                painter.setPen(QColor(col))
                if base:
                    painter.setFont(self._get_math_font(base_pt, italic=True))
                    painter.drawText(int(cur_x), int(y_center + 4), base)
                    cur_x += QFontMetrics(painter.font()).horizontalAdvance(base)
                painter.setFont(self._get_math_font(max(base_pt - 3, 8)))
                painter.drawText(int(cur_x), int(y_center - 2), sup)
                cur_x += QFontMetrics(painter.font()).horizontalAdvance(sup) + 2

            elif kind == 'exp':
                _, exp_txt, col = tok
                painter.setPen(QColor(col))
                painter.setFont(self._get_math_font(base_pt, italic=True))
                fm_e = QFontMetrics(painter.font())
                painter.drawText(int(cur_x), int(y_center + 4), "e")
                cur_x += fm_e.horizontalAdvance("e") + 1
                painter.setFont(self._get_math_font(max(base_pt - 3, 9)))
                painter.drawText(int(cur_x), int(y_center - 2), exp_txt)
                cur_x += QFontMetrics(painter.font()).horizontalAdvance(exp_txt) + 3

            elif kind == 'frac':
                _, num_toks, den_toks, col = tok
                cur_x = self._draw_frac(painter, cur_x, y_center, num_toks, den_toks, col, base_pt)

        return cur_x


# ==============================================================================
# 3. ESQUEMÁTICO INTERACTIVO CON FÍSICA Y PARTÍCULAS POR TRAMO DE CABLE
# ==============================================================================
class CircuitCanvas(QWidget):
    def __init__(self, engine: TransientEngine, on_switch_changed_callback):
        super().__init__()
        self.engine = engine
        self.on_switch_changed = on_switch_changed_callback
        self.setMinimumHeight(245)
        self.setMouseTracking(True)

        # Fases de partículas independientes por cada conductor físico
        self.PARTICLE_STEP = 28.0
        self.ph_source = 0.0
        self.ph_r1 = 0.0
        self.ph_r2 = 0.0
        self.ph_load_top = 0.0
        self.ph_L = 0.0
        self.ph_C = 0.0
        self.ph_load_bot = 0.0
        self.ph_gnd_return = 0.0

        self.sw1_rect = QRectF()
        self.sw2_rect = QRectF()

    @staticmethod
    def _calc_step(curr):
        if abs(curr) < 1e-5:
            return 0.0
        s = np.clip(curr * 28.0, -12.0, 12.0)
        if abs(s) < 1.0:
            s = 1.0 if curr > 0 else -1.0
        return s

    def advance_particles(self):
        s = self.PARTICLE_STEP
        self.ph_source = (self.ph_source + self._calc_step(self.engine.i_source)) % s
        self.ph_r1 = (self.ph_r1 + self._calc_step(self.engine.i_source)) % s
        self.ph_r2 = (self.ph_r2 + self._calc_step(self.engine.iR2)) % s
        self.ph_load_top = (self.ph_load_top + self._calc_step(self.engine.i_load)) % s
        self.ph_L = (self.ph_L + self._calc_step(self.engine.iL)) % s
        self.ph_C = (self.ph_C + self._calc_step(self.engine.iC)) % s
        self.ph_load_bot = (self.ph_load_bot + self._calc_step(self.engine.i_load)) % s
        self.ph_gnd_return = (self.ph_gnd_return + self._calc_step(self.engine.i_source)) % s
        self.update()

    def mouseMoveEvent(self, event):
        pos = event.pos()
        if self.sw1_rect.contains(pos) or self.sw2_rect.contains(pos):
            self.setCursor(Qt.PointingHandCursor)
        else:
            self.setCursor(Qt.ArrowCursor)

    def mousePressEvent(self, event):
        pos = event.pos()
        if self.sw1_rect.contains(pos):
            self.engine.sw1_state = (self.engine.sw1_state + 1) % 4
            self.on_switch_changed("S1 (Clic Esquemático)")
        elif self.sw2_rect.contains(pos):
            self.engine.sw2_closed = not self.engine.sw2_closed
            self.on_switch_changed("S2 (Clic Esquemático)")

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)

        w = self.width()
        h = self.height()

        painter.fillRect(0, 0, w, h, QColor("#0D111A"))
        grid_pen = QPen(QColor("#182130"), 1)
        painter.setPen(grid_pen)
        for x in range(0, w, 24):
            painter.drawLine(x, 0, x, h)
        for y in range(0, h, 24):
            painter.drawLine(0, y, w, y)

        top_y = int(h * 0.23)
        bot_y = int(h * 0.79)
        mid_y = int((top_y + bot_y) / 2)

        left_x = int(w * 0.09)
        sw1_x = int(w * 0.23)
        r1_x = int(w * 0.38)
        r2_x = int(w * 0.52)
        load_x1 = int(w * 0.70)
        load_x2 = int(w * 0.88)

        wire_pen = QPen(QColor("#475569"), 3)
        painter.setPen(wire_pen)

        # Cable inferior común de masa
        painter.drawLine(left_x, bot_y, load_x2, bot_y)
        self._draw_ground(painter, r2_x, bot_y)

        # =====================================================================
        # 1. BLOQUE DE FUENTES Y LLAVE S1 (4 ESTADOS)
        # =====================================================================
        self.sw1_rect = QRectF(sw1_x - 42, top_y - 30, 90, 60)
        painter.setPen(QPen(QColor("#0284C7"), 1, Qt.DashLine))
        painter.setBrush(QBrush(QColor(14, 116, 144, 35)))
        painter.drawRoundedRect(self.sw1_rect, 6, 6)

        painter.setPen(wire_pen)
        painter.drawLine(left_x, top_y + 18, left_x, bot_y)

        st = self.engine.sw1_state
        if st == TransientEngine.SW_DIRECT:
            self._draw_source(painter, left_x, mid_y, f"+{self.engine.Vs:.0f}V", inverted=False)
            painter.setPen(QPen(QColor("#00E5FF"), 3))
            painter.drawLine(left_x, top_y + 18, left_x, top_y)
            painter.drawLine(left_x, top_y, sw1_x + 22, top_y)
            s1_txt = f"S1: Directa (+{self.engine.Vs:.0f}V)"
            s1_col = QColor("#00E5FF")

        elif st == TransientEngine.SW_INVERT:
            self._draw_source(painter, left_x, mid_y, f"-{self.engine.Vs:.0f}V", inverted=True)
            painter.setPen(QPen(QColor("#FB7185"), 3))
            painter.drawLine(left_x, top_y + 18, left_x, top_y)
            painter.drawLine(left_x, top_y, sw1_x + 22, top_y)
            s1_txt = f"S1: Invertida (-{self.engine.Vs:.0f}V)"
            s1_col = QColor("#FB7185")

        elif st == TransientEngine.SW_SHORT:
            self._draw_source(painter, left_x, mid_y, "0V Corto", inverted=False, dimmed=True)
            painter.setPen(QPen(QColor("#34D399"), 3))
            painter.drawLine(sw1_x - 12, bot_y, sw1_x - 12, top_y)
            painter.drawLine(sw1_x - 12, top_y, sw1_x + 22, top_y)
            s1_txt = "S1: Corto (0V Homog.)"
            s1_col = QColor("#34D399")

        else:
            self._draw_source(painter, left_x, mid_y, "Abierta", inverted=False, dimmed=True)
            painter.setPen(QPen(QColor("#FBBF24"), 3))
            painter.drawLine(left_x, top_y, sw1_x - 10, top_y)
            painter.drawLine(sw1_x - 10, top_y, sw1_x + 16, top_y - 20)
            s1_txt = "S1: Abierta (Homog. Pura)"
            s1_col = QColor("#FBBF24")

        painter.setFont(QFont("Segoe UI", 9, QFont.Bold))
        painter.setPen(s1_col)
        painter.drawText(int(sw1_x - 42), int(top_y - 34), f"[Clic] {s1_txt}")

        # =====================================================================
        # 2. RESISTENCIA PRINCIPAL R1
        # =====================================================================
        painter.setPen(wire_pen)
        painter.drawLine(sw1_x + 22, top_y, r1_x - 26, top_y)
        p1 = (self.engine.i_source ** 2) * self.engine.R1
        self._draw_resistor_h(painter, r1_x, top_y, f"R1={self.engine.R1:.0f}Ω", p1)
        painter.drawLine(r1_x + 26, top_y, r2_x, top_y)

        # =====================================================================
        # 3. RAMA TRANSVERSAL S2 + R2
        # =====================================================================
        self.sw2_rect = QRectF(r2_x - 30, top_y + 4, 60, 42)
        painter.setPen(QPen(QColor("#D97706"), 1, Qt.DashLine))
        painter.setBrush(QBrush(QColor(120, 53, 15, 45)))
        painter.drawRoundedRect(self.sw2_rect, 6, 6)

        painter.setPen(wire_pen)
        painter.drawLine(r2_x, top_y, r2_x, top_y + 12)
        painter.setPen(QPen(QColor("#FBBF24"), 3))
        if self.engine.sw2_closed:
            painter.drawLine(r2_x, top_y + 12, r2_x, top_y + 38)
            s2_txt = "S2: Cerrada"
        else:
            painter.drawLine(r2_x, top_y + 12, r2_x + 16, top_y + 32)
            s2_txt = "S2: Abierta"

        painter.setFont(QFont("Segoe UI", 8, QFont.Bold))
        painter.setPen(QColor("#FBBF24"))
        painter.drawText(r2_x - 30, top_y - 8, f"[Clic] {s2_txt}")

        painter.setPen(wire_pen)
        painter.drawLine(r2_x, top_y + 38, r2_x, bot_y)
        p2 = (self.engine.iR2 ** 2) * self.engine.R2
        self._draw_resistor_v(painter, r2_x, int((top_y + 38 + bot_y) / 2), f"R2={self.engine.R2:.0f}Ω", p2)

        # =====================================================================
        # 4. RED REACTIVA SEGÚN TOPOLOGÍA
        # =====================================================================
        painter.setPen(wire_pen)
        if self.engine.topology == 'RLC_SERIE':
            painter.drawLine(r2_x, top_y, load_x1 - 25, top_y)
            self._draw_inductor_h(painter, load_x1, top_y, f"L={self.engine.format_eng(self.engine.L, 'H')}", self.engine.iL)
            painter.drawLine(load_x1 + 25, top_y, load_x2, top_y)
            painter.drawLine(load_x2, top_y, load_x2, bot_y)
            self._draw_capacitor_v(painter, load_x2, mid_y, f"C={self.engine.format_eng(self.engine.C, 'F')}", self.engine.vC)

        elif self.engine.topology == 'RC_SERIE':
            painter.drawLine(r2_x, top_y, load_x2, top_y)
            painter.drawLine(load_x2, top_y, load_x2, bot_y)
            self._draw_capacitor_v(painter, load_x2, mid_y, f"C={self.engine.format_eng(self.engine.C, 'F')}", self.engine.vC)

        elif self.engine.topology == 'RL_SERIE':
            painter.drawLine(r2_x, top_y, load_x2, top_y)
            painter.drawLine(load_x2, top_y, load_x2, bot_y)
            self._draw_inductor_v(painter, load_x2, mid_y, f"L={self.engine.format_eng(self.engine.L, 'H')}", self.engine.iL)

        elif self.engine.topology == 'RLC_PARALELO':
            painter.drawLine(r2_x, top_y, load_x2, top_y)
            # Rama Inductor L
            painter.drawLine(load_x1, top_y, load_x1, bot_y)
            self._draw_inductor_v(painter, load_x1, mid_y, f"L={self.engine.format_eng(self.engine.L, 'H')}", self.engine.iL)
            # Rama Capacitor C
            painter.drawLine(load_x2, top_y, load_x2, bot_y)
            self._draw_capacitor_v(painter, load_x2, mid_y, f"C={self.engine.format_eng(self.engine.C, 'F')}", self.engine.vC)

        # =====================================================================
        # 5. ANIMACIÓN CONTINUA DE PARTÍCULAS POR CADA CONDUCTOR FÍSICO
        # =====================================================================
        # A) Rama de Fuente y R1
        if st in (TransientEngine.SW_DIRECT, TransientEngine.SW_INVERT):
            self._draw_particles(painter, left_x, bot_y, left_x, top_y, self.ph_source, self.engine.i_source)
            self._draw_particles(painter, left_x, top_y, sw1_x + 22, top_y, self.ph_source, self.engine.i_source)
            self._draw_particles(painter, sw1_x + 22, top_y, r2_x, top_y, self.ph_r1, self.engine.i_source)
            self._draw_particles(painter, r2_x, bot_y, left_x, bot_y, self.ph_gnd_return, self.engine.i_source)

        elif st == TransientEngine.SW_SHORT:
            self._draw_particles(painter, sw1_x - 12, bot_y, sw1_x - 12, top_y, self.ph_source, self.engine.i_source)
            self._draw_particles(painter, sw1_x - 12, top_y, sw1_x + 22, top_y, self.ph_source, self.engine.i_source)
            self._draw_particles(painter, sw1_x + 22, top_y, r2_x, top_y, self.ph_r1, self.engine.i_source)
            self._draw_particles(painter, r2_x, bot_y, sw1_x - 12, bot_y, self.ph_gnd_return, self.engine.i_source)

        # B) Rama transversal R2
        if self.engine.sw2_closed:
            self._draw_particles(painter, r2_x, top_y, r2_x, bot_y, self.ph_r2, self.engine.iR2)

        # C) Red Reactiva según Topología
        if self.engine.topology == 'RLC_SERIE':
            self._draw_particles(painter, r2_x, top_y, load_x2, top_y, self.ph_L, self.engine.iL)
            self._draw_particles(painter, load_x2, top_y, load_x2, bot_y, self.ph_L, self.engine.iL)
            self._draw_particles(painter, load_x2, bot_y, r2_x, bot_y, self.ph_load_bot, self.engine.iL)

        elif self.engine.topology == 'RC_SERIE':
            self._draw_particles(painter, r2_x, top_y, load_x2, top_y, self.ph_C, self.engine.iC)
            self._draw_particles(painter, load_x2, top_y, load_x2, bot_y, self.ph_C, self.engine.iC)
            self._draw_particles(painter, load_x2, bot_y, r2_x, bot_y, self.ph_load_bot, self.engine.iC)

        elif self.engine.topology == 'RL_SERIE':
            self._draw_particles(painter, r2_x, top_y, load_x2, top_y, self.ph_L, self.engine.iL)
            self._draw_particles(painter, load_x2, top_y, load_x2, bot_y, self.ph_L, self.engine.iL)
            self._draw_particles(painter, load_x2, bot_y, r2_x, bot_y, self.ph_load_bot, self.engine.iL)

        elif self.engine.topology == 'RLC_PARALELO':
            self._draw_particles(painter, r2_x, top_y, load_x1, top_y, self.ph_load_top, self.engine.i_load)
            self._draw_particles(painter, load_x1, top_y, load_x1, bot_y, self.ph_L, self.engine.iL)
            self._draw_particles(painter, load_x1, top_y, load_x2, top_y, self.ph_C, self.engine.iC)
            self._draw_particles(painter, load_x2, top_y, load_x2, bot_y, self.ph_C, self.engine.iC)
            self._draw_particles(painter, load_x2, bot_y, load_x1, bot_y, self.ph_C, self.engine.iC)
            self._draw_particles(painter, load_x1, bot_y, r2_x, bot_y, self.ph_load_bot, self.engine.i_load)

        # =====================================================================
        # 6. LECTURAS INSTANTÁNEAS EN EL ESQUEMÁTICO
        # =====================================================================
        painter.setFont(QFont("Consolas", 10, QFont.Bold))
        fm = painter.fontMetrics()

        v_tag = "v_C" if self.engine.topology != 'RL_SERIE' else "v_L"
        txt_v = f"{v_tag}(t) = {self.engine.vC:+.2f} V"

        i_tag = "i_L" if self.engine.topology != 'RC_SERIE' else "i_C"
        i_val = self.engine.iL if self.engine.topology != 'RC_SERIE' else self.engine.iC
        txt_i = f"{i_tag}(t) = {i_val*1e3:+.2f} mA"

        txt_src = f"i_fuente = {self.engine.i_source*1e3:+.2f} mA   │   R_eq = {self.engine.format_eng(self.engine.r_eq, 'Ω')}"

        w_src = fm.horizontalAdvance(txt_src)
        w_i = fm.horizontalAdvance(txt_i)
        w_v = fm.horizontalAdvance(txt_v)

        painter.setPen(QColor("#38BDF8"))
        painter.drawText(left_x, bot_y + 24, txt_src)

        x_i = max(left_x + w_src + 35, int(w * 0.53))
        painter.setPen(QColor("#FFEA00"))
        painter.drawText(x_i, bot_y + 24, txt_i)

        x_v = max(x_i + w_i + 30, int(w * 0.77))
        painter.setPen(QColor("#00E676"))
        painter.drawText(x_v, bot_y + 24, txt_v)

    def _draw_particles(self, painter, x1, y1, x2, y2, phase, current_val):
        if abs(current_val) < 1e-5:
            return
        dist = math.hypot(x2 - x1, y2 - y1)
        if dist < 2:
            return
        ux = (x2 - x1) / dist
        uy = (y2 - y1) / dist

        painter.setPen(Qt.NoPen)
        painter.setBrush(QBrush(QColor("#FDE047")))
        step = self.PARTICLE_STEP
        d = phase % step
        while d < dist:
            px = x1 + ux * d
            py = y1 + uy * d
            painter.drawEllipse(QPointF(px, py), 3.8, 3.8)
            d += step

    def _draw_ground(self, painter, x, y):
        painter.setPen(QPen(QColor("#94A3B8"), 2))
        painter.drawLine(x, y, x, y + 8)
        painter.drawLine(x - 12, y + 8, x + 12, y + 8)
        painter.drawLine(x - 8, y + 12, x + 8, y + 12)
        painter.drawLine(x - 4, y + 16, x + 4, y + 16)

    def _draw_source(self, painter, x, y, label, inverted=False, dimmed=False):
        col = QColor("#475569") if dimmed else (QColor("#FB7185") if inverted else QColor("#00E5FF"))
        painter.setPen(QPen(col, 2))
        painter.setBrush(QBrush(QColor("#0D111A")))
        painter.drawEllipse(QPointF(x, y), 18, 18)
        painter.setFont(QFont("Consolas", 10, QFont.Bold))
        top_sign = "−" if inverted else "+"
        bot_sign = "+" if inverted else "−"
        painter.drawText(x - 4, y - 3, top_sign)
        painter.drawText(x - 4, y + 12, bot_sign)
        painter.setFont(QFont("Consolas", 8, QFont.Bold))
        painter.drawText(x - 78, y + 5, label)

    def _draw_resistor_h(self, painter, x, y, label, power=0.0):
        heat_ratio = np.clip(power / 2.0, 0.0, 1.0)
        core_color = QColor(
            int(30 + 190 * heat_ratio),
            int(38 + 50 * heat_ratio * 0.4),
            int(54 * (1.0 - heat_ratio * 0.7))
        )
        border_color = QColor("#FF9100") if heat_ratio < 0.3 else QColor("#FF3D00")

        painter.setPen(QPen(border_color, 2))
        painter.setBrush(QBrush(core_color))
        painter.drawRect(x - 26, y - 9, 52, 18)
        painter.setFont(QFont("Consolas", 8, QFont.Bold))
        painter.setPen(QColor("#CBD5E1"))
        painter.drawText(x - 24, y - 14, label)
        if power > 1e-4:
            painter.setPen(QColor("#FF7043"))
            painter.setFont(QFont("Segoe UI", 7))
            p_str = f"P={power*1e3:.0f}mW" if power < 1.0 else f"P={power:.2f}W"
            painter.drawText(x - 22, y + 22, p_str)

    def _draw_resistor_v(self, painter, x, y, label, power=0.0):
        heat_ratio = np.clip(power / 2.0, 0.0, 1.0)
        core_color = QColor(
            int(30 + 190 * heat_ratio),
            int(38 + 50 * heat_ratio * 0.4),
            int(54 * (1.0 - heat_ratio * 0.7))
        )
        border_color = QColor("#FF9100") if heat_ratio < 0.3 else QColor("#FF3D00")

        painter.setPen(QPen(border_color, 2))
        painter.setBrush(QBrush(core_color))
        painter.drawRect(x - 9, y - 22, 18, 44)
        painter.setFont(QFont("Consolas", 8, QFont.Bold))
        painter.setPen(QColor("#CBD5E1"))
        painter.drawText(x + 14, y + 4, label)
        if power > 1e-4:
            painter.setPen(QColor("#FF7043"))
            painter.setFont(QFont("Segoe UI", 7))
            p_str = f"P={power*1e3:.0f}mW" if power < 1.0 else f"P={power:.2f}W"
            painter.drawText(x + 14, y + 17, p_str)

    def _draw_inductor_h(self, painter, x, y, label, current=0.0):
        mag_strength = np.clip(abs(current) / 0.1, 0.0, 1.0)
        if mag_strength > 0.05:
            glow_color = QColor(255, 234, 0, int(60 * mag_strength))
            painter.setPen(QPen(glow_color, 2, Qt.DashLine))
            painter.setBrush(Qt.NoBrush)
            painter.drawEllipse(QPointF(x, y), 36, 18)

        painter.setPen(QPen(QColor("#FFEA00"), 2))
        painter.setBrush(QBrush(QColor("#0D111A")))
        for i in range(-2, 2):
            painter.drawEllipse(QPointF(x + i * 12 + 6, y), 7, 9)
        painter.setFont(QFont("Consolas", 8, QFont.Bold))
        painter.setPen(QColor("#CBD5E1"))
        painter.drawText(x - 25, y - 14, label)

    def _draw_inductor_v(self, painter, x, y, label, current=0.0):
        mag_strength = np.clip(abs(current) / 0.1, 0.0, 1.0)
        if mag_strength > 0.05:
            glow_color = QColor(255, 234, 0, int(60 * mag_strength))
            painter.setPen(QPen(glow_color, 2, Qt.DashLine))
            painter.setBrush(Qt.NoBrush)
            painter.drawEllipse(QPointF(x, y), 18, 36)

        painter.setPen(QPen(QColor("#FFEA00"), 2))
        painter.setBrush(QBrush(QColor("#0D111A")))
        for i in range(-2, 2):
            painter.drawEllipse(QPointF(x, y + i * 12 + 6), 9, 7)
        painter.setFont(QFont("Consolas", 8, QFont.Bold))
        painter.setPen(QColor("#CBD5E1"))
        painter.drawText(x - 80, y + 4, label)

    def _draw_capacitor_v(self, painter, x, y, label, voltage=0.0):
        painter.setPen(QPen(QColor("#0D111A"), 6))
        painter.drawLine(x, y - 8, x, y + 8)

        charge_ratio = np.clip(abs(voltage) / max(abs(self.engine.Vs), 1.0), 0.0, 1.0)
        glow_alpha = int(70 * charge_ratio)
        if glow_alpha > 5:
            field_brush = QColor(0, 230, 118, glow_alpha) if voltage >= 0 else QColor(56, 189, 248, glow_alpha)
            painter.fillRect(x - 16, y - 7, 32, 14, field_brush)

        plate_color = QColor("#00E676") if voltage >= 0 else QColor("#38BDF8")
        painter.setPen(QPen(plate_color, 4))
        painter.drawLine(x - 16, y - 7, x + 16, y - 7)
        painter.drawLine(x - 16, y + 7, x + 16, y + 7)

        painter.setFont(QFont("Consolas", 8, QFont.Bold))
        if voltage > 1.0:
            painter.setPen(QColor("#00E676"))
            painter.drawText(x + 20, y - 3, "+ +")
            painter.drawText(x + 20, y + 13, "− −")
        elif voltage < -1.0:
            painter.setPen(QColor("#38BDF8"))
            painter.drawText(x + 20, y - 3, "− −")
            painter.drawText(x + 20, y + 13, "+ +")

        painter.setPen(QColor("#CBD5E1"))
        painter.drawText(x - 82, y - 12, label)


# ==============================================================================
# 4. VENTANA PRINCIPAL DE LA APLICACIÓN
# ==============================================================================
class RLCSimulatorApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Simulador Interactivo RLC (v2.0) — Ecuaciones Diferenciales, Régimen y Respuesta Homogénea")
        self.resize(1400, 920)
        self.setStyleSheet("""
            QMainWindow, QWidget { background-color: #0A0D15; color: #E2E8F0; font-family: 'Segoe UI', sans-serif; }
            QGroupBox { border: 1px solid #233044; border-radius: 6px; margin-top: 8px; font-weight: bold; color: #38BDF8; padding-top: 8px; font-size: 11.5px; }
            QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 5px; }
            QPushButton { background-color: #1E293B; border: 1px solid #475569; border-radius: 5px; padding: 6px 10px; font-weight: bold; color: #F8FAFC; text-align: left; font-size: 11px; }
            QPushButton:hover { background-color: #334155; border-color: #38BDF8; }
            QPushButton:pressed { background-color: #0F172A; }
            QDoubleSpinBox, QComboBox { background-color: #151D2C; border: 1px solid #3B4B64; border-radius: 4px; padding: 4px 6px; color: #FFFFFF; font-weight: bold; font-size: 11px; }
        """)

        self.engine = TransientEngine()
        self.paused = False
        self.tau_rate = 1.2

        # Buffer del osciloscopio
        self.max_pts = 600
        self.vC_buf = np.zeros(self.max_pts)
        self.iL_buf = np.zeros(self.max_pts)
        self.vth_buf = np.zeros(self.max_pts)

        # Marcadores verticales de conmutación
        self.transient_markers = []

        self._build_ui()
        self.update_source_buttons_ui()
        self.update_s2_button_ui()
        self.add_transient_marker("t₀ (Inicio)")
        self.refresh_all_views()

        # Timer de simulación (50 FPS)
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.on_timer_tick)
        self.timer.start(20)

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QHBoxLayout(central)
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setSpacing(10)

        # --- PANEL IZQUIERDO ---
        left_panel = QWidget()
        left_panel.setFixedWidth(380)
        left_layout = QVBoxLayout(left_panel)
        left_layout.setSpacing(8)
        left_layout.setContentsMargins(0, 0, 0, 0)

        # 1. Topología
        grp_topo = QGroupBox("1. Topología del Circuito")
        topo_layout = QVBoxLayout(grp_topo)
        self.combo_topo = QComboBox()
        self.combo_topo.addItems([
            "RLC Serie (2° Orden)",
            "RLC Paralelo (2° Orden)",
            "RC Serie (1° Orden)",
            "RL Serie (1° Orden)"
        ])
        self.combo_topo.currentIndexChanged.connect(self.on_topology_changed)
        topo_layout.addWidget(self.combo_topo)
        left_layout.addWidget(grp_topo)

        # 2. Componentes
        grp_comp = QGroupBox("2. Valores de Componentes")
        grid_comp = QGridLayout(grp_comp)
        grid_comp.setHorizontalSpacing(8)
        grid_comp.setVerticalSpacing(5)

        self.spin_Vs = self._make_spin(1.0, 200.0, 12.0, " V", 1.0)
        self.spin_R1 = self._make_spin(1.0, 10000.0, 100.0, " Ω", 10.0)
        self.spin_R2 = self._make_spin(1.0, 10000.0, 220.0, " Ω", 10.0)
        self.spin_L = self._make_spin(0.1, 5000.0, 50.0, " mH", 5.0)
        self.spin_C = self._make_spin(0.1, 5000.0, 10.0, " µF", 1.0)

        grid_comp.addWidget(QLabel("Tensión Fuente (|Vs|):"), 0, 0)
        grid_comp.addWidget(self.spin_Vs, 0, 1)
        grid_comp.addWidget(QLabel("Resistencia R1 (Serie):"), 1, 0)
        grid_comp.addWidget(self.spin_R1, 1, 1)
        grid_comp.addWidget(QLabel("Resistencia R2 (Paralelo):"), 2, 0)
        grid_comp.addWidget(self.spin_R2, 2, 1)
        grid_comp.addWidget(QLabel("Inductancia (L):"), 3, 0)
        grid_comp.addWidget(self.spin_L, 3, 1)
        grid_comp.addWidget(QLabel("Capacitancia (C):"), 4, 0)
        grid_comp.addWidget(self.spin_C, 4, 1)
        left_layout.addWidget(grp_comp)

        # 3. Control de Régimen y Conmutación (Forzada vs Homogénea)
        grp_trans = QGroupBox("3. Régimen y Disparo de Transitorios")
        trans_layout = QVBoxLayout(grp_trans)
        trans_layout.setSpacing(5)

        self.lbl_state_badge = QLabel("EN TRANSITORIO")
        self.lbl_state_badge.setAlignment(Qt.AlignCenter)
        self.lbl_state_badge.setStyleSheet("background-color: #B45309; color: #FEF3C7; font-weight: bold; padding: 6px; border-radius: 4px; font-size: 11px;")
        trans_layout.addWidget(self.lbl_state_badge)

        btn_steady = QPushButton("⚡ Forzar Régimen Permanente (t → ∞)")
        btn_steady.setStyleSheet("background-color: #065F46; border-color: #10B981; color: #A7F3D0; text-align: center;")
        btn_steady.clicked.connect(self.on_force_steady)
        trans_layout.addWidget(btn_steady)

        self.btn_direct = QPushButton("🔌 1. Fuente Directa (+Vs) — [EDO Forzada]")
        self.btn_direct.clicked.connect(lambda checked=False: self.set_s1(TransientEngine.SW_DIRECT))
        trans_layout.addWidget(self.btn_direct)

        self.btn_invert = QPushButton("🔄 2. Invertir Fuente (-Vs) — [EDO Forzada]")
        self.btn_invert.clicked.connect(lambda checked=False: self.set_s1(TransientEngine.SW_INVERT))
        trans_layout.addWidget(self.btn_invert)

        self.btn_short = QPushButton("📉 3. Fuente en Corto (0V) — [EDO Homogénea]")
        self.btn_short.clicked.connect(lambda checked=False: self.set_s1(TransientEngine.SW_SHORT))
        trans_layout.addWidget(self.btn_short)

        self.btn_open = QPushButton("🔓 4. Fuente Abierta (∞Ω) — [Solo Homogénea]")
        self.btn_open.clicked.connect(lambda checked=False: self.set_s1(TransientEngine.SW_OPEN))
        trans_layout.addWidget(self.btn_open)

        self.btn_s2 = QPushButton("🔀 Llave S2: CERRADA (R2 activa)")
        self.btn_s2.clicked.connect(self.toggle_s2)
        trans_layout.addWidget(self.btn_s2)

        btn_zero = QPushButton("⚡ Descargar Variables a Cero (vC=0, iL=0)")
        btn_zero.setStyleSheet("background-color: #334155; border-color: #64748B; color: #CBD5E1; text-align: center;")
        btn_zero.clicked.connect(self.on_reset_zero)
        trans_layout.addWidget(btn_zero)

        left_layout.addWidget(grp_trans)

        # 4. Escalas y Velocidad
        grp_time = QGroupBox("4. Velocidad Temporal y Escala Fija")
        time_layout = QGridLayout(grp_time)

        self.combo_speed = QComboBox()
        self.combo_speed.addItems([
            "Auto-Sync Lenta (0.5 τ / seg)",
            "Auto-Sync Normal (1.2 τ / seg)",
            "Auto-Sync Rápida (2.5 τ / seg)",
            "Auto-Sync Ultra (5.0 τ / seg)"
        ])
        self.combo_speed.setCurrentIndex(1)
        self.combo_speed.currentIndexChanged.connect(self.on_speed_changed)

        self.spin_vmax = self._make_spin(2.0, 500.0, 20.0, " V", 5.0, trigger_transient=False)
        self.spin_vmax.valueChanged.connect(self.update_plot_axes)

        btn_pause = QPushButton("⏸ Pausar / Reanudar Animación")
        btn_pause.setStyleSheet("text-align: center;")
        btn_pause.clicked.connect(self.toggle_pause)

        time_layout.addWidget(QLabel("Velocidad Relativa a τ:"), 0, 0)
        time_layout.addWidget(self.combo_speed, 0, 1)
        time_layout.addWidget(QLabel("Escala Vertical Fija (±V):"), 1, 0)
        time_layout.addWidget(self.spin_vmax, 1, 1)
        time_layout.addWidget(btn_pause, 2, 0, 1, 2)
        left_layout.addWidget(grp_time)
        left_layout.addStretch()

        # --- PANEL DERECHO ---
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.setSpacing(8)

        # 1. Esquemático
        self.canvas = CircuitCanvas(self.engine, self.on_canvas_switch_clicked)
        right_layout.addWidget(self.canvas, stretch=3)

        # 2. Pizarra de Ecuaciones Diferenciales en formato grande
        self.math_widget = MathEquationWidget(self.engine)
        right_layout.addWidget(self.math_widget, stretch=3)

        # 3. Osciloscopio con ejes fijos bloqueados
        self.plot_widget = pg.PlotWidget(title="Osciloscopio en Tiempo Real (Ejes Bloqueados)")
        self.plot_widget.setBackground("#080B12")
        self.plot_widget.showGrid(x=True, y=True, alpha=0.25)

        # Bloqueo estricto de mouse y rueda
        self.plot_widget.setMouseEnabled(x=False, y=False)
        self.plot_widget.hideButtons()
        self.plot_widget.setMenuEnabled(False)

        self.plot_widget.setLabel('left', 'Tensión [V]  /  Corriente Escalada (i · R_eq) [V]')
        self.plot_widget.setLabel('bottom', 'Ventana Temporal Sincronizada con τ')
        self.plot_widget.setYRange(-20, 20, padding=0)
        self.plot_widget.setXRange(0, self.max_pts, padding=0)

        self.plot_widget.addLegend(offset=(12, 12))
        self.curve_vth = self.plot_widget.plot(pen=pg.mkPen('#FF9100', width=1.5, style=Qt.DashLine), name="Excitación V_th [V]")
        self.curve_vC = self.plot_widget.plot(pen=pg.mkPen('#00E676', width=2.6), name="Tensión v_C(t) / v_L(t) [V]")
        self.curve_iL = self.plot_widget.plot(pen=pg.mkPen('#FFEA00', width=2.2), name="Corriente i(t) · R_eq [V]")

        right_layout.addWidget(self.plot_widget, stretch=4)

        main_layout.addWidget(left_panel)
        main_layout.addWidget(right_panel, stretch=1)

    def _make_spin(self, min_v, max_v, val, suffix, step, trigger_transient=True):
        sp = QDoubleSpinBox()
        sp.setRange(min_v, max_v)
        sp.setDecimals(2)
        sp.setValue(val)
        sp.setSuffix(suffix)
        sp.setSingleStep(step)
        if trigger_transient:
            sp.valueChanged.connect(self.on_params_changed)
        return sp

    def update_source_buttons_ui(self):
        st = self.engine.sw1_state
        btns = [self.btn_direct, self.btn_invert, self.btn_short, self.btn_open]
        active_styles = [
            "background-color: #0284C7; border: 2px solid #38BDF8; color: #FFFFFF; font-weight: bold; text-align: left;",
            "background-color: #BE123C; border: 2px solid #FB7185; color: #FFFFFF; font-weight: bold; text-align: left;",
            "background-color: #047857; border: 2px solid #34D399; color: #FFFFFF; font-weight: bold; text-align: left;",
            "background-color: #B45309; border: 2px solid #FBBF24; color: #FFFFFF; font-weight: bold; text-align: left;",
        ]
        default_style = "background-color: #1E293B; border: 1px solid #475569; color: #94A3B8; font-weight: bold; text-align: left;"
        for i, b in enumerate(btns):
            if i == st:
                b.setStyleSheet(active_styles[i])
            else:
                b.setStyleSheet(default_style)

    def update_s2_button_ui(self):
        if self.engine.sw2_closed:
            self.btn_s2.setText("🔀 Llave S2: CERRADA (R2 activa)")
            self.btn_s2.setStyleSheet("background-color: #78350F; border: 1.5px solid #F59E0B; color: #FEF3C7; font-weight: bold; text-align: left;")
        else:
            self.btn_s2.setText("🔀 Llave S2: ABIERTA (R2 desconectada)")
            self.btn_s2.setStyleSheet("background-color: #1E293B; border: 1px solid #475569; color: #94A3B8; font-weight: bold; text-align: left;")

    def add_transient_marker(self, label="t₀"):
        line = pg.InfiniteLine(
            pos=self.max_pts - 1,
            angle=90,
            pen=pg.mkPen('#00E5FF', width=1.5, style=Qt.DashLine),
            label=label,
            labelOpts={'position': 0.92, 'color': (0, 229, 255), 'movable': False, 'fill': (15, 23, 42, 200)}
        )
        self.plot_widget.addItem(line)
        self.transient_markers.append({'line': line, 'x': float(self.max_pts - 1)})

    def on_topology_changed(self, idx):
        mapping = ['RLC_SERIE', 'RLC_PARALELO', 'RC_SERIE', 'RL_SERIE']
        self.engine.topology = mapping[idx]
        self.engine.recalculate_transient()
        self.add_transient_marker("t₀ (Topología)")
        self.refresh_all_views()

    def on_params_changed(self):
        self.engine.Vs = self.spin_Vs.value()
        self.engine.R1 = self.spin_R1.value()
        self.engine.R2 = self.spin_R2.value()
        self.engine.L = self.spin_L.value() * 1e-3
        self.engine.C = self.spin_C.value() * 1e-6
        self.engine.recalculate_transient()
        self.add_transient_marker("t₀ (Parám.)")
        self.refresh_all_views()

    def set_s1(self, state):
        self.engine.sw1_state = state
        self.engine.recalculate_transient()
        names = ["+Vs", "-Vs", "0V Corto", "Abierta"]
        self.add_transient_marker(f"t₀ ({names[state]})")
        self.update_source_buttons_ui()
        self.refresh_all_views()

    def toggle_s2(self):
        self.engine.sw2_closed = not self.engine.sw2_closed
        self.engine.recalculate_transient()
        self.add_transient_marker("t₀ (S2)")
        self.update_s2_button_ui()
        self.refresh_all_views()

    def on_canvas_switch_clicked(self, tag):
        self.engine.recalculate_transient()
        self.add_transient_marker(f"t₀ ({tag})")
        self.update_source_buttons_ui()
        self.update_s2_button_ui()
        self.refresh_all_views()

    def on_force_steady(self):
        self.engine.force_steady_state()
        self.add_transient_marker("t→∞ (Forzado)")
        self.refresh_all_views()

    def on_reset_zero(self):
        self.engine.vC = 0.0
        self.engine.iL = 0.0
        self.engine.recalculate_transient()
        self.add_transient_marker("t₀ (Reset 0)")
        self.refresh_all_views()

    def on_speed_changed(self, idx):
        rates = [0.5, 1.2, 2.5, 5.0]
        self.tau_rate = rates[idx]

    def update_plot_axes(self):
        vmax = self.spin_vmax.value()
        self.plot_widget.setYRange(-vmax, vmax, padding=0)

    def toggle_pause(self):
        self.paused = not self.paused

    def refresh_all_views(self):
        self.math_widget.update()
        self.canvas.update()

    def on_timer_tick(self):
        if self.paused:
            return

        dt = 0.020 * self.tau_rate * max(self.engine.tau, 1e-9)
        self.engine.step(dt)
        self.canvas.advance_particles()

        taus_elapsed = self.engine.t_local / max(self.engine.tau, 1e-9)
        if self.engine.is_in_steady_state():
            self.lbl_state_badge.setText(f"RÉGIMEN PERMANENTE ESTABLECIDO (t' = {taus_elapsed:.1f} τ > 5τ)")
            self.lbl_state_badge.setStyleSheet("background-color: #065F46; color: #A7F3D0; font-weight: bold; padding: 6px; border-radius: 4px; font-size: 11px;")
        else:
            self.lbl_state_badge.setText(f"EN TRANSITORIO — Evolución: {taus_elapsed:.2f} τ  (Régimen a las 5τ)")
            self.lbl_state_badge.setStyleSheet("background-color: #B45309; color: #FEF3C7; font-weight: bold; padding: 6px; border-radius: 4px; font-size: 11px;")

        # Desplazar marcadores verticales del osciloscopio
        alive_markers = []
        for m in self.transient_markers:
            m['x'] -= 1.0
            if m['x'] >= 0:
                m['line'].setValue(m['x'])
                alive_markers.append(m)
            else:
                self.plot_widget.removeItem(m['line'])
        self.transient_markers = alive_markers

        # Desplazar muestras del osciloscopio
        self.vC_buf[:-1] = self.vC_buf[1:]
        self.iL_buf[:-1] = self.iL_buf[1:]
        self.vth_buf[:-1] = self.vth_buf[1:]

        # Escala de corriente multiplicada por R_eq (unidades de Volt)
        r_scale = self.engine.R1 if math.isinf(self.engine.r_eq) else self.engine.r_eq
        i_plot = self.engine.iC if self.engine.topology == 'RC_SERIE' else self.engine.iL

        self.vC_buf[-1] = self.engine.vC
        self.iL_buf[-1] = i_plot * r_scale
        self.vth_buf[-1] = self.engine.v_th

        x_axis = np.arange(self.max_pts)
        self.curve_vth.setData(x_axis, self.vth_buf)
        self.curve_vC.setData(x_axis, self.vC_buf)
        self.curve_iL.setData(x_axis, self.iL_buf)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = RLCSimulatorApp()
    win.show()
    sys.exit(app.exec_())
