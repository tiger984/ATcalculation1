class Line:
    def __init__(self,
                 name="",                # 导线名
                 type_name="",           # 导线型号名
                 resistance=0,           # 导线电阻
                 radius=1,               # 导线计算半径
                 equivalent_radius=1,    # 导线等效半径
                 rho=1,                  # 导线电导率
                 mu_r=1,                 # 相对磁导率
                 coordinate_x=0,         # 导线坐标 x
                 coordinate_y=0):        # 导线坐标 y
        self.name = name
        self.type_name = type_name
        self.resistance = resistance
        self.radius = radius
        self.rho = rho
        self.mu_r = mu_r
        self.equivalent_radius = equivalent_radius
        self.coordinate_x = coordinate_x
        self.coordinate_y = coordinate_y

    def set_parameter(self,
                      resistance=0,
                      radius=1,
                      equivalent_radius=1,
                      rho=1,
                      mu_r=1,
                      coordinate_x=0,
                      coordinate_y=0
                      ):
        self.resistance = resistance
        self.radius = radius
        self.rho = rho
        self.mu_r = mu_r
        self.equivalent_radius = equivalent_radius
        self.coordinate_x = coordinate_x
        self.coordinate_y = coordinate_y


class ATLines:
    lines_type = "AT System"

    def __init__(self,
                 name=""
                 ):
        self.name = name
        self.lines = []
        self.lines.append(Line(name="cw1"))
        self.lines.append(Line(name="mw1"))
        self.lines.append(Line(name="pf1"))
        self.lines.append(Line(name="ra1"))
        self.lines.append(Line(name="ra2"))
        self.lines.append(Line(name="pw1"))
        self.lines.append(Line(name="e1"))
        self.lines.append(Line(name="cw2"))
        self.lines.append(Line(name="mw2"))
        self.lines.append(Line(name="pf2"))
        self.lines.append(Line(name="ra3"))
        self.lines.append(Line(name="ra4"))
        self.lines.append(Line(name="pw2"))
        self.lines.append(Line(name="e2"))

   #  def __setattr__(self, file_name=""):
        # for line in lines:
        #     pass
    # def set_accodinater(self):
     #        pass


if __name__ == '__main__':
    ATlines = ATLines(name="AT")
    n = len(ATlines.lines)
    print(ATlines.name)
    for i in range(n):
        print(ATlines.lines[i].name, ATlines.lines[i].coordinate_x)


class TractionTransformer:               # 牵引变压器类
    def __init__(self,
                 name="",
                 alias_name="",          # 牵引变压器型号
                 zs=1+1j,                # 系统内阻抗
                 location=0,             # 安装位置（公里)

                 ):
        self.name = name
        self.alias_name = alias_name
        self.zs = zs
        self.location = location


class AutoTransformer:                    # 自耦变压器类
    def __init__(self,
                 name="",
                 alias_name="",           # 变压器型号
                 zs=1+1j,                 # 漏抗（欧）
                 location=0             # 安装位置(公里）
                 ):
        self.name = name
        self.alias_name = alias_name
        self.zs = zs
        self.location = location


class Locomotive:                        # 机车类
    def __init__(self,
                 name="",                #
                 load=0,                 # 负荷电流（安 ,A）
                 harmonic={},            # 谐波含量 （fn, %）
                 location=0,             # 运行位置 （公里，km）
                 at_upline=1             # 位于上行线=1；位于下行线=0
                 ):
        self.name = name
        self.load = load
        self.harmonic = harmonic
        self.location = location


class CrossConnection:
    """

    """
    def __init__(self,
                 name="",
                 ):
        self.name = name
        self.to_all = []           # 上下行全并联，( km）
        self.pw1_ra1 = []          # 保护线pw1连钢轨ra1，( km）
        self.e1_ra1 = []           # 综合地线e1连接钢轨ra1，( km）
        self.e1_g = {}             # 综合地线e1连接大地g ，{ km, 欧}
        self.pw2_ra3 = []          # 保护线pw2连钢轨ra3，( km）
        self. e2_ra3 = []          # 综合地线e2连接钢轨ra3，( km）
        self.e2_g = {}             # 综合地线e2连接大地g，{ km, 欧}
        self.ra1_g = {}            # 钢轨ra1连接大地g
        self.ra3_g = {}            # 钢轨ra3连接大地g


class Topology:
    def __init__(self,
                 name="",
                 alias_name="",
                 section_length=25.6,
                 earth_rou=10**6,
                 traction_transformer=[],
                 lines=[],
                 auto_transformers=[],
                 cross_connections=[]
                 ):
        self.name = name
        self.alias_name = alias_name
        self.section_length = section_length
        self.earth_rou = earth_rou
        self.traction_transformer = traction_transformer
        self.lines = lines
        self.auto_transformers = auto_transformers
        self.cross_connections = cross_connections

    def set_topology(self,
                     db_file_name="default.db"    # 系统拓扑结构 sqlite数据库文件名
                     ):
        pass


# 设定测试topology原始数据
cross_connections = CrossConnection()
cross_connections.to_all = [0, 15, 25.6]
cross_connections.ra1_g = {0: 0.5, 15: 0.9, 15.6: 0.9}
cross_connections.ra2_g = {0: 0.5, 15: 0.9, 15.6: 0.9}
cross_connections.e1_g = {0: 0.5, 15: 0.9, 15.6: 0.9}
cross_connections.e2_g = {0: 0.5, 15: 0.9, 15.6: 0.9}
cross_connections.

topology = Topology(name="测试供电臂")
topology.cross_connections = cross_connections

