import matplotlib.pyplot as plt
import numpy as np

bcc_energies=[-55.2651986628596, -55.70054835374419, -55.80087740016136, -55.78256362811685, -55.72558295286919, -55.66319836996898, -55.60512698519512]
bcc_volumes=[35.9281, 52.7214, 74.0697, 100.5198, 132.6183, 170.9117, 215.9467]

hcp_energies= [-55.141581951558564, -55.705768166476716, -55.783408147749995, -55.76843343425499, -55.70773864062073]
hcp_volumes=np.array([32.72975, 50.81005, 74.55925, 104.7504, 142.1565])


plt.figure()
plt.plot(bcc_volumes, bcc_energies, label = 'bcc')
plt.plot(hcp_volumes, hcp_energies, label = 'hcp')
plt.plot([62.5, 62.5], [-55.85, -55.15], linestyle = '--', c = 'red', linewidth = 3)
plt.xlabel('Volume (a.u.^3)')
plt.ylabel('Energy (eV)')
plt.legend()
plt.xlim(25, 225)
plt.ylim(-55.85, -55.15)
plt.savefig('1B_bcc_hcp.png', dpi = 100)
plt.show()