// Enable toggle
function toggleDarkMode() {
  document.documentElement.classList.toggle('dark');
  localStorage.setItem('theme', document.documentElement.classList.contains('dark') ? 'dark' : 'light');
}

// Apply saved theme
if (localStorage.getItem('theme') === 'dark') {
  document.documentElement.classList.add('dark');
}