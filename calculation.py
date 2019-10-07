# calculation.py
#
import numpy as np
from scipy import linalg
import topology as tp
import conductors_calculation as cc



class ChainNetwork:

    def __init__(self,
                 name="",
                 h=1,    # 基波、谐波次数。h = 1,基波；h >1 ,谐波
                 delta_length=0.5,    # 基本分段长度（km）
                 m=6,    # 归并后的导线数
                 topology=tp.Topology(name="AT_System"),   # 供电网络拓扑结构数据
                 ):

        self.name = name
        self.topology = topology
        self.h = h
        self.m = m
        self.section_length = self.topology.section_length    # 供电臂总长度
        self.delta_length = delta_length    # 基本分段长度
        self.n = int(np.round(self.section_length / self.delta_length) + 1)  # 分割面数，即网络节点数
        self.distances = self.__set_distances()                # 分割点位置距离
        self.segment_lengths = self.__set_segment_lengths()    # 各分割段长度

        self.z = np.zeros((self.m, self.m, self.n-1), np.complex128)    # 网络纵向阻抗矩阵
        self.y = np.zeros((self.m, self.m, self.n), np.complex128)      # 网络横向导纳矩阵
        self.Un = np.zeros((self.m, self.n), np.complex128)     # 节点电压向量
        self.In = np.zeros((self.m, self.n-1), np.complex128)   # 纵向导线电流向量
        self.Gn = np.zeros((self.m, self.n), np.complex128)     # 节点注入电流向量


    def reset(self, topology, h, m, delta_length):
        self.topology = topology
        self.h = h
        self.m = m
        self.delta_length = delta_length
        self.section_length = self.topology.section_length
        self.name = self.topology.name
        self.n = int(np.round(self.section_length / self.delta_length) + 1)  # 分割面数，即网络节点
        self.distances = self.__set_distances()  # 分割点位置距离
        self.segment_lengths = self.__set_segment_lengths()  # 各分割段长度

    def __set_distances(self):
        distances = np.arange(0, self.n, 1) * self.delta_length
        distances[self.n - 1] = self.section_length
        return distances

    def __set_segment_lengths(self):
        segment_lengths = np.zeros(self.n-1, np.float64)
        segment_lengths[:n-1] = self.distances[1:self.n] - self.distances[0:self.n-1]    # 分割段长度 （km）
        # print("segment_num = ", self.segment_lengths.shape)
        return segment_lengths

    def set_z_y(self):
        unit_z = self.__calc_unit_z()
        unit_yc = self.__calc_unit_yc()

        for i in range(self.n-1):
            self.z[:, :, i] = unit_z*self.segment_lengths[i]
        for i in range(1, self.n-1):
            self.y[:, :, i] = unit_yc * (self.segment_lengths[i-1]/2 + self.segment_lengths[i]/2)
            self.y[:, :, 0] = unit_yc * (self.segment_lengths[0]/2)
            self.y[:, :, self.n-1] = unit_yc * (self.segment_lengths[self.n-2]/2)

        print("y=", self.y[0, 0, 0])
        print(self.z.shape)


    def construct_M(self):
        n = self.n
        m = self.m
        Dk = np.zeros((m, m, n), np.complex128)
        Mk = np.zeros((m, m, n), np.complex128)
        M = np.zeros((m*n, m*n), np.complex128)
        for k in range(1, n):   # 生产Dk元素矩阵
            Dk[:, :, k] = -linalg.inv(self.z[:, :, k-1])

        for k in range(1, n-1):    # 生成Mk元素矩阵
            Mk[:, :, k] = self.y[:, :, k] + linalg.inv(self.z[:, :, k]) + linalg.inv(self.z[:, :, k - 1])

        Mk[:, :, 0] = self.y[:, :, 0] + linalg.inv(self.z[:, :, 0])
        Mk[:, :, n-1] = self.y[:, :, n-1] + linalg.inv(self.z[:, :, n-2])

        for i in range(n):
            if i == 0:
                self.M[0:m, 0:m] = Mk[:, :, 0]
            else:
                k = i*m  # 分块起点号
                M[k:k+m, k-m:k] = Dk[:, :, i]
                M[k - m: k, k:k + m] = Dk[:, :, i]
                M[k:k+m, k:k+m] = Mk[:, :, i]
        return M

    def solution(self):
       # self.set_z_y()
       # self.add_y()
        # #self.add_G()
        M = self.construct_M()
        W = linalg.inv(self.M)
        G = self.Gn.reshape(self.m*self.n)
        U = np.dot(W, G)
        Un = U.reshape(self.m, self.n)
        for k in range(self.n-1):
            z = self.z[:, :, k]
            z = linalg.inv(z)
            u = Un[:, k+1]-Un[:, k]
            In[:, k] = np.dot(z, u)

    def add_y(self):
        cross_connections = self.topology.cross_connections
        y0 = 10**6    # 连接导线的电导率
        if self.m == 6:
            for x in cross_connections.to_all:        # 全并联连接线
                k = int(np.round(x/self.delta_length))
                self.y[:, :, k] = self.y[:, :, k] + y6_to_all*y0

            for x in cross_connections.ra1_g:
                k = int(np.round(x/self.delta_length))
                if k < self.n:
                    rg = cross_connections.ra1_g.get(x)
                    self.y[:, :, k] = self.y[:, :, k] + y6_ra1_g * (1/rg)

            for x in cross_connections.ra3_g:
                k = int(np.round(x / self.delta_length))
                if k < self.n:
                    rg = cross_connections.ra3_g.get(x)
                    self.y[:, :, k] = self.y[:, :, k] + y6_ra3_g * (1 / rg)

            for x in cross_connections.e1_g:
                k = int(np.round(x / self.delta_length))
                if k < self.n:
                    rg = cross_connections.e1_g.get(x)
                    self.y[:, :, k] = self.y[:, :, k] + y6_e1_g * (1 / rg)

            for x in cross_connections.e2_g:
                k = int(np.round(x / self.delta_length))
                if k < self.n:
                    rg = cross_connections.e2_g.get(x)
                    self.y[:, :, k] = self.y[:, :, k] + y6_e2_g * (1 / rg)

    def __get_line_parameter(self):    # 获取导线参数
        c_xy = []
        c_resistance = []
        c_mu_r = []
        c_radius = []
        c_equivalent_radius = []
        c_rho = []
        earth_rou = self.topology.earth_rou
        for line in self.topology.lines:
            c_xy.append([line.coordinater_x, line.coordinater_y])
            c_resistance.append(line.resistance)
            c_mu_r.append(line.mu_r)
            c_radius.append(line.radius)
            c_equivalent_radius.append(line.equivalent_radius)
            c_rho.append(line.rho)
        return c_xy, earth_rou, c_resistance, c_mu_r, c_radius, c_equivalent_radius, c_rho

    def __get_line_parameter_test(self):    # 设定测试导线参数
        c_xy = cc.conductors_coordinator
        c_resistance = cc.Rd
        c_mu_r = cc.mu_r
        c_radius = cc.conductors_calc_radius
        c_equivalent_radius = cc.conductors_equivalent_radius
        c_rho = cc.rho
        earth_rou = self.topology.earth_rou
        return c_xy, earth_rou, c_resistance, c_mu_r, c_radius, c_equivalent_radius, c_rho

    def __calc_unit_z(self):
        """
        计算平行导线原始单位阻抗矩阵
        :return:     阻抗矩阵
        """
        c_xy, earth_rou, c_resistance, c_mu_r, c_radius, c_equivalent_radius, c_rho = self.__get_line_parameter_test()
        f = 50*self.h     # 频率
        z = cc.calc_Zf(f, c_xy, c_radius, c_resistance, c_rho, c_mu_r, earth_rou)
        if self.m == 6:
            z = cc.merge_z(z, 0, 1)
            z = cc.merge_z(z, 2, 3)
            z = cc.merge_z(z, 2, 3)
            z = cc.merge_z(z, 2, 3)
            z = cc.merge_z(z, 3, 4)
            z = cc.merge_z(z, 5, 6)
            z = cc.merge_z(z, 5, 6)
            z = cc.merge_z(z, 5, 6)
        unit_z = z
        return unit_z

    def __calc_unit_yc(self):
        """
        计算导线的单位电容导纳矩阵
        :return:
        """
        c_xy, earth_rou, c_resistance, c_mu_r, c_radius, c_equivalent_radius, c_rho = self.__get_line_parameter_test()
        p = cc.calc_potential_coefficient(c_xy, c_radius)
        if self.m == 6:
            p = cc.merge_potential_coefficient(p, 0, 1)    # 合并成6根导线
            p = cc.merge_potential_coefficient(p, 2, 3)
            p = cc.merge_potential_coefficient(p, 2, 3)
            p = cc.merge_potential_coefficient(p, 2, 3)
            p = cc.merge_potential_coefficient(p, 3, 4)
            p = cc.merge_potential_coefficient(p, 5, 6)
            p = cc.merge_potential_coefficient(p, 5, 6)
            p = cc.merge_potential_coefficient(p, 5, 6)
        b = cc.calc_B(p)    # 计算电容
        f = 50*self.h
        unit_yc = -1j*2*np.pi*f*b   # 计算电容导纳
        return unit_yc


# 定义基本导纳关联矩阵 (m=6)

AT_matrix_6 = cc.set_connction_matrix_AT(6, 0, 1, 2)  # AT导纳关联矩阵
y6_to_all = cc.set_connection_matrix(6, 0, 3) + cc.set_connection_matrix(6, 1, 4) \
                 + cc.set_connection_matrix(6, 2, 5)  # 上下行全并联导纳连接矩阵
y6_e1_g = cc.set_connction_matrix_g(6, 1)  # 综合地线e1连接大地g ，导纳关联矩阵
y6_e2_g = cc.set_connction_matrix_g(6, 4)  # 综合地线e2连接大地g ，导纳关联矩
y6_ra1_g = cc.set_connction_matrix_g(6, 1)  # 钢轨ra1连接大地g, 导纳关联矩阵
y6_ra3_g = cc.set_connction_matrix_g(6, 4)  # 钢轨ra3连接大地g, 导纳关联矩阵


if __name__ == '__main__':

    chen = ChainNetwork(h=1, m=6, delta_length=0.1)
    chen.set_z_y()
    chen.construct_M()

    print(chen.M.shape)

    # print(chen.AT_matrix_6)
    # print(chen.y6_to_all)
    # print(chen.y6_e1_g)