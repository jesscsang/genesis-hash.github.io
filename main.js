const btn = document.getElementById('btn');

btn.addEventListener('click', function handleClick() {
    if (btn.textContent == 'Click Me') {
        btn.textContent = 'Stop';
        btn.style.background = 'red';
    }
    else {
        btn.textContent = 'Click Me';
        btn.style.background = 'green';
    }
});
